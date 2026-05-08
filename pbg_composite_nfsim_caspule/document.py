"""Build the Composite document from the wiring map.

The three configurable inputs are the three input files the user
edits to drive a run:

    caspule_input_file: path to a LAMMPS / CASPULE ``.in`` script
    nfsim_model_file:   path to a BNGL ``.bngl`` file
    detector_config:    path to the detector YAML

Plus per-process intervals. CASPULE's interval should be much smaller
than NFSim's: CASPULE runs the spatial dynamics that the detector
watches, while NFSim integrates a separate non-spatial pool that only
needs to react when the detector injects matter.
"""

from pbg_composite_nfsim_caspule.detector import (
    DetectorConfig,
    load_detector_config,
    _parse_detector_config,
)
from pbg_composite_nfsim_caspule.wiring import WIRING, SCALAR_EMIT_TYPES


def _resolve_detector_config(
    detector_config: str | None,
    detector_rules: list | None,
) -> DetectorConfig:
    if detector_config:
        return load_detector_config(detector_config)
    if detector_rules:
        return _parse_detector_config({'rules': detector_rules})
    raise ValueError(
        'build_document needs detector_config (path) or detector_rules (list)')


def build_document(
    *,
    caspule_input_file: str = '',
    caspule_input_script: str = '',
    nfsim_model_file: str,
    detector_config: str = '',
    detector_rules: list | None = None,
    caspule_interval: float = 0.05,
    nfsim_interval: float = 0.5,
    nfsim_n_steps: int = 50,
    decoupled: bool = False,
    emit_scalars: bool = True,
):
    """Return a Composite document wiring caspule + detector + nfsim.

    Args:
        caspule_input_file: path to a LAMMPS .in script. Mutually
            exclusive with ``caspule_input_script``.
        caspule_input_script: inline LAMMPS / CASPULE script.
        nfsim_model_file: path to a BNGL model.
        detector_config: path to a YAML rule file. Mutually exclusive
            with ``detector_rules``.
        detector_rules: inline list of rule dicts.
        caspule_interval: CASPULE step in LAMMPS time units. Default
            0.05 lj-time-units (~10 timesteps at dt=0.005).
        nfsim_interval: NFSim step in seconds. Default 0.5.
        nfsim_n_steps: NFSim ``n_steps`` per BNGL ``simulate(...)`` call.
        decoupled: if True, route the detector's outputs to dead-end
            stores so neither CASPULE nor NFSim receive them. Used by
            the demo as a baseline that proves coupling matters.
        emit_scalars: attach the RAM emitter for the SCALAR_EMIT_TYPES
            ports plus species. Default True.

    Returns:
        A document dict suitable for ``Composite({'state': doc}, core=...)``.
    """
    if not caspule_input_file and not caspule_input_script:
        raise ValueError(
            'build_document requires caspule_input_file or caspule_input_script')

    detector_cfg = _resolve_detector_config(detector_config, detector_rules)
    obs_names = detector_cfg.observable_names

    # Override wiring for the decoupled baseline: route the detector's
    # OUTPUTS to dead-end stores. Consumers (caspule.atoms_to_remove,
    # nfsim.observables) stay on the normal stores — they just never
    # see anything because nobody else writes to those stores in this
    # composition. Net effect: CASPULE keeps all atoms, NFSim's species
    # pool stays empty, but every process still has well-typed inputs.
    if decoupled:
        atoms_to_remove_out = ['stores', '_decoupled', 'atoms_to_remove']
        species_add_path    = ['stores', '_decoupled', 'species_additions']
    else:
        atoms_to_remove_out = WIRING[('caspule', 'atoms_to_remove')]
        species_add_path    = WIRING[('detector', 'molecule_additions')]
    atoms_to_remove_in = WIRING[('caspule', 'atoms_to_remove')]
    species_in_path    = WIRING[('nfsim', 'observables_in')]

    caspule_outputs = {
        port: WIRING[('caspule', port)]
        for port in (
            'temperature', 'potential_energy', 'kinetic_energy',
            'total_energy', 'pressure', 'volume', 'box_dimensions',
            'num_atoms', 'positions', 'velocities', 'atom_types',
            'num_bonds', 'bonds', 'bonds_by_type', 'bond_energy',
            'formed_bonds', 'broken_bonds',
            'num_clusters', 'largest_cluster', 'cluster_sizes',
        )
    }

    doc = {
        # --- CASPULE: bond-aware MD ---
        'caspule': {
            '_type': 'process',
            'address': 'local:CASPULEProcess',
            'config': {
                'input_file': caspule_input_file,
                'input_script': caspule_input_script,
            },
            'interval': caspule_interval,
            'inputs': {
                'atoms_to_remove': atoms_to_remove_in,
            },
            'outputs': caspule_outputs,
        },

        # --- ObservableDetector: bridges spatial -> non-spatial ---
        'detector': {
            '_type': 'step',
            'address': 'local:ObservableDetector',
            'config': {
                'config_file': detector_config,
                'rules': detector_rules or [],
            },
            'inputs': {
                'positions':   WIRING[('detector', 'positions')],
                'bonds':       WIRING[('detector', 'bonds')],
                'atom_types':  WIRING[('detector', 'atom_types')],
                'num_atoms':   WIRING[('detector', 'num_atoms')],
            },
            'outputs': {
                'atoms_to_remove':    atoms_to_remove_out,
                'molecule_additions': species_add_path,
            },
        },

        # --- NFSim: rule-based non-spatial kinetics ---
        'nfsim': {
            '_type': 'process',
            'address': 'local:NFSimProcess',
            'config': {
                'model_file': nfsim_model_file,
                'n_steps': nfsim_n_steps,
            },
            'interval': nfsim_interval,
            'inputs':  {'observables': species_in_path},
            'outputs': {'observables': WIRING[('nfsim', 'observables_out')]},
        },

        # --- shared stores ---
        'stores': {
            'positions': [],
            'bonds': [],
            'atom_types': [],
            'num_atoms': 0,
            'atoms_to_remove': [],
            'species': {name: 0.0 for name in obs_names},
            '_decoupled': {
                'atoms_to_remove': [],
                'species_additions': {name: 0.0 for name in obs_names},
            },
        },
    }

    if emit_scalars:
        emit_schema = dict(SCALAR_EMIT_TYPES)
        emit_schema['time'] = 'float'
        emit_schema['species'] = 'map[float]'
        emitter_inputs = {
            p: WIRING[('caspule', p)] for p in SCALAR_EMIT_TYPES
        }
        emitter_inputs['time'] = ['global_time']
        emitter_inputs['species'] = ['stores', 'species']
        doc['emitter'] = {
            '_type': 'step',
            'address': 'local:ram-emitter',
            'config': {'emit': emit_schema},
            'inputs': emitter_inputs,
        }

    return doc
