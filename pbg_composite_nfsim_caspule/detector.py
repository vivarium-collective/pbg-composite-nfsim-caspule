"""ObservableDetector Step.

Reads CASPULE's spatial state (positions, bonds, atom_types) and
converts user-defined detection rules into two outputs:

  - ``atoms_to_remove``: a list of CASPULE global atom IDs to delete
    on the next CASPULE step.
  - ``molecule_additions``: a ``map[float]`` of NFSim observable name ->
    delta count to inject into NFSim's non-spatial pool.

A rule says "find clusters / atoms in CASPULE matching this pattern,
each match counts as one molecule of <nfsim_observable>, and (optionally)
remove the matched atoms from CASPULE". The detector is configured by a
YAML file (the input-file standard for this Step) listing rules.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import Optional

import yaml
from process_bigraph import Step


# -- rule data classes --------------------------------------------------------

@dataclass
class DetectorRule:
    """A single detection rule.

    Exactly one of ``bond_type_cluster`` or ``atom_type`` is the primary
    selector; the other ``match_*`` fields constrain it further.

    Attributes:
        name: NFSim observable name to credit each match to.
        bond_type_cluster: if set, find connected components of bonds
            of this LAMMPS bond-type integer; each component is a match.
        atom_type: if set, find unbonded atoms of this LAMMPS atom-type;
            each atom is a match.
        min_size: cluster size threshold (only applies to bond_type_cluster).
        max_size: cluster size upper bound (only applies to bond_type_cluster).
        max_bonds: when ``atom_type`` is set, only count atoms whose bond
            count is at-most this. ``0`` selects free monomers.
        remove_from_caspule: whether matched atoms are removed from
            CASPULE after detection.
    """

    name: str
    bond_type_cluster: Optional[int] = None
    atom_type: Optional[int] = None
    min_size: int = 1
    max_size: Optional[int] = None
    max_bonds: Optional[int] = None
    remove_from_caspule: bool = True

    def __post_init__(self):
        if (self.bond_type_cluster is None) == (self.atom_type is None):
            raise ValueError(
                f'rule {self.name!r}: set exactly one of '
                'bond_type_cluster or atom_type')
        if self.min_size < 1:
            raise ValueError(f'rule {self.name!r}: min_size must be >= 1')


@dataclass
class DetectorConfig:
    """A list of DetectorRule entries.

    The observable_names property returns the NFSim observables this
    detector will write to — the document builder uses it to seed the
    species store with zeros.
    """

    rules: list[DetectorRule] = field(default_factory=list)

    @property
    def observable_names(self) -> list[str]:
        return [r.name for r in self.rules]


def load_detector_config(path: str) -> DetectorConfig:
    """Load a DetectorConfig from a YAML file.

    Schema:
        rules:
          - name: str
            (bond_type_cluster: int | atom_type: int)
            min_size: int = 1
            max_size: int | null = null
            max_bonds: int | null = null
            remove_from_caspule: bool = true
    """
    with open(path) as f:
        data = yaml.safe_load(f) or {}
    return _parse_detector_config(data)


def _parse_detector_config(data: dict) -> DetectorConfig:
    raw_rules = data.get('rules', [])
    rules = [DetectorRule(**r) for r in raw_rules]
    return DetectorConfig(rules=rules)


# -- Step implementation ------------------------------------------------------

class ObservableDetector(Step):
    """Step that turns CASPULE spatial state into NFSim species counts.

    Configured by either ``config_file`` (path to YAML) or ``rules``
    (inline list of dicts). On every CASPULE write, the detector
    computes the connected-component graph from the bond list, applies
    each rule, and emits both the matched atom IDs and the per-rule
    counts as deltas.

    The ``atoms_to_remove`` output uses ``overwrite[list]`` because the
    detector is the sole authority each tick on what to delete, and a
    delta-list would compose nonsensically.
    """

    config_schema = {
        'config_file': {'_type': 'string', '_default': ''},
        'rules': {'_type': 'list', '_default': []},
    }

    def __init__(self, config=None, core=None):
        super().__init__(config=config, core=core)
        if self.config['config_file']:
            self._cfg = load_detector_config(self.config['config_file'])
        elif self.config['rules']:
            self._cfg = _parse_detector_config({'rules': self.config['rules']})
        else:
            raise ValueError(
                'ObservableDetector requires either config_file or rules')

    @property
    def observable_names(self) -> list[str]:
        return self._cfg.observable_names

    def inputs(self):
        return {
            'positions': 'list',
            'bonds': 'list',
            'atom_types': 'list',
            'num_atoms': 'integer',
        }

    def outputs(self):
        return {
            'atoms_to_remove': 'overwrite[list]',
            'molecule_additions': 'map[float]',
        }

    def update(self, state):
        bonds = state.get('bonds') or []
        atom_types = state.get('atom_types') or []
        num_atoms = int(state.get('num_atoms') or 0)

        atoms_to_remove: set[int] = set()
        additions: dict[str, float] = {r.name: 0.0 for r in self._cfg.rules}

        # Per-atom bond counts (independent of bond type). Used by atom_type
        # rules that gate on max_bonds.
        bond_count_by_atom: dict[int, int] = {}
        for b in bonds:
            _btype, a1, a2 = int(b[0]), int(b[1]), int(b[2])
            bond_count_by_atom[a1] = bond_count_by_atom.get(a1, 0) + 1
            bond_count_by_atom[a2] = bond_count_by_atom.get(a2, 0) + 1

        for rule in self._cfg.rules:
            if rule.bond_type_cluster is not None:
                clusters = _components_of_bond_type(
                    bonds, rule.bond_type_cluster)
                for atoms in clusters:
                    if len(atoms) < rule.min_size:
                        continue
                    if rule.max_size is not None and len(atoms) > rule.max_size:
                        continue
                    additions[rule.name] += 1.0
                    if rule.remove_from_caspule:
                        atoms_to_remove.update(atoms)
            else:
                # atom_type rule
                target_type = rule.atom_type
                # atom_types is a 1D list indexed by 0-based local position.
                # CASPULE emits global atom IDs implicitly via order: index
                # i corresponds to atom ID i+1 in the LAMMPS sense (since
                # gather_bonds always uses global IDs). We recover the
                # mapping by enumerating types in order.
                for i, t in enumerate(atom_types):
                    if int(t) != target_type:
                        continue
                    atom_id = i + 1
                    if rule.max_bonds is not None:
                        if bond_count_by_atom.get(atom_id, 0) > rule.max_bonds:
                            continue
                    additions[rule.name] += 1.0
                    if rule.remove_from_caspule:
                        atoms_to_remove.add(atom_id)

        # CASPULE expects an authoritative list each tick (overwrite[list]).
        # Sorting keeps test output deterministic.
        return {
            'atoms_to_remove': sorted(atoms_to_remove),
            'molecule_additions': additions,
        }


# -- helpers ------------------------------------------------------------------

def _components_of_bond_type(bonds: list, bond_type: int) -> list[list[int]]:
    """Return connected components of the subgraph induced by bonds of
    type ``bond_type``. Each component is a sorted list of atom IDs."""
    edges = [(int(a), int(b)) for (t, a, b) in bonds if int(t) == bond_type]
    if not edges:
        return []

    # Union-find on the atoms touched by these edges.
    parent: dict[int, int] = {}

    def find(x: int) -> int:
        parent.setdefault(x, x)
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> None:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in edges:
        union(a, b)

    groups: dict[int, list[int]] = {}
    for atom in list(parent.keys()):
        r = find(atom)
        groups.setdefault(r, []).append(atom)
    return [sorted(g) for g in groups.values()]
