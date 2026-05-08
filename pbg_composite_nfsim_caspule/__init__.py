"""Composite of CASPULE (bond-aware MD) + NFSim (rule-based kinetics)
coupled through an observable detector that converts spatial bond
clusters in CASPULE into species counts in NFSim's non-spatial pool.

Public surface:
    build_core()        - allocate a process-bigraph core with every class registered
    build_document()    - build the Composite document from the wiring map
    ObservableDetector  - the detector Step
    DetectorConfig      - parsed YAML rule list

Each simulator is configured by a standard input file:
    CASPULE  <- LAMMPS .in script
    NFSim    <- BNGL .bngl file
    Detector <- YAML rules file
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

__all__ = [
    'ObservableDetector',
    'DetectorConfig',
    'DetectorRule',
    'load_detector_config',
    'build_core',
    'build_document',
    'WIRING',
]
