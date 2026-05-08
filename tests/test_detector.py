"""Unit tests for ObservableDetector — pure-function tests on update()
without spinning up a full Composite. The detector is the only piece
of glue between CASPULE and NFSim, so we test it in isolation first."""

import os
import textwrap

import pytest

from pbg_composite_nfsim_caspule.core import build_core
from pbg_composite_nfsim_caspule.detector import (
    DetectorRule,
    ObservableDetector,
    load_detector_config,
)


@pytest.fixture
def core():
    return build_core()


def _make_detector(core, rules):
    return ObservableDetector(config={'rules': rules}, core=core)


def test_dimer_rule_emits_one_per_bond_type_cluster(core):
    """Two disjoint type-1 bonds → two Dimer matches → +2 in NFSim."""
    rules = [{'name': 'Dimer', 'bond_type_cluster': 1,
              'min_size': 2, 'max_size': 2}]
    det = _make_detector(core, rules)
    state = {
        'positions':  [[0, 0, 0]] * 4,
        'bonds':      [[1, 1, 2], [1, 3, 4]],
        'atom_types': [1, 1, 1, 1],
        'num_atoms':  4,
    }
    out = det.update(state)
    assert out['molecule_additions'] == {'Dimer': 2.0}
    assert out['atoms_to_remove'] == [1, 2, 3, 4]


def test_dimer_rule_skips_too_large_clusters(core):
    """A 3-atom chain (bonds 1-2 and 2-3) should not register as a Dimer."""
    rules = [{'name': 'Dimer', 'bond_type_cluster': 1,
              'min_size': 2, 'max_size': 2}]
    det = _make_detector(core, rules)
    state = {
        'positions':  [[0, 0, 0]] * 3,
        'bonds':      [[1, 1, 2], [1, 2, 3]],
        'atom_types': [1, 1, 1],
        'num_atoms':  3,
    }
    out = det.update(state)
    assert out['molecule_additions'] == {'Dimer': 0.0}
    assert out['atoms_to_remove'] == []


def test_atom_type_rule_with_max_bonds_finds_free_monomers(core):
    """`max_bonds: 0` selects atoms that are not in any bond."""
    rules = [{'name': 'FreeMonomer', 'atom_type': 1, 'max_bonds': 0,
              'remove_from_caspule': False}]
    det = _make_detector(core, rules)
    state = {
        'positions':  [[0, 0, 0]] * 4,
        'bonds':      [[1, 2, 3]],
        'atom_types': [1, 1, 1, 1],
        'num_atoms':  4,
    }
    out = det.update(state)
    # atoms 1 and 4 are unbonded -> 2 free monomers
    assert out['molecule_additions'] == {'FreeMonomer': 2.0}
    assert out['atoms_to_remove'] == []  # remove_from_caspule=False


def test_remove_from_caspule_false_does_not_remove(core):
    rules = [{'name': 'Dimer', 'bond_type_cluster': 1,
              'min_size': 2, 'max_size': 2,
              'remove_from_caspule': False}]
    det = _make_detector(core, rules)
    state = {
        'positions':  [[0, 0, 0]] * 2,
        'bonds':      [[1, 1, 2]],
        'atom_types': [1, 1],
        'num_atoms':  2,
    }
    out = det.update(state)
    assert out['molecule_additions'] == {'Dimer': 1.0}
    assert out['atoms_to_remove'] == []


def test_yaml_config_round_trips(core, tmp_path):
    cfg_path = tmp_path / 'rules.yaml'
    cfg_path.write_text(textwrap.dedent("""
        rules:
          - name: Dimer
            bond_type_cluster: 1
            min_size: 2
            max_size: 2
          - name: FreeMonomer
            atom_type: 1
            max_bonds: 0
            remove_from_caspule: false
    """))
    cfg = load_detector_config(str(cfg_path))
    assert [r.name for r in cfg.rules] == ['Dimer', 'FreeMonomer']
    assert cfg.rules[0].bond_type_cluster == 1
    assert cfg.rules[1].remove_from_caspule is False


def test_invalid_rule_neither_selector_raises(core):
    with pytest.raises(ValueError, match='exactly one of'):
        DetectorRule(name='X')


def test_invalid_rule_both_selectors_raises(core):
    with pytest.raises(ValueError, match='exactly one of'):
        DetectorRule(name='X', bond_type_cluster=1, atom_type=1)
