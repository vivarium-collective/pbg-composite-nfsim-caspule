"""Composite of CASPULE (bond-aware MD) + NFSim (rule-based kinetics)
coupled through an observable detector that converts spatial bond
clusters in CASPULE into species counts in NFSim's non-spatial pool.

The fastest way to run a simulation is to point the bundled CLI at a
PBG JSON document plus optional overrides for the three input files::

    python -m pbg_composite_nfsim_caspule examples/composite_dimer.pbg.json

Or in Python::

    from pbg_composite_nfsim_caspule import run
    run('examples/composite_dimer.pbg.json', total_time=2.0)

Each simulator is configured by its own input file standard:

    CASPULE  <- LAMMPS .in script
    NFSim    <- BNGL .bngl file
    Detector <- YAML rules file

The PBG JSON references those by relative path, resolved against the
JSON file's directory.
"""

from pbg_composite_nfsim_caspule.detector import (
    ObservableDetector,
    DetectorConfig,
    DetectorRule,
    load_detector_config,
)
from pbg_composite_nfsim_caspule.core import build_core
from pbg_composite_nfsim_caspule.document import build_document
from pbg_composite_nfsim_caspule.wiring import WIRING
from pbg_composite_nfsim_caspule.run import (
    run,
    load_document,
    resolve_input_paths,
)

__all__ = [
    'ObservableDetector',
    'DetectorConfig',
    'DetectorRule',
    'load_detector_config',
    'build_core',
    'build_document',
    'WIRING',
    'run',
    'load_document',
    'resolve_input_paths',
]
