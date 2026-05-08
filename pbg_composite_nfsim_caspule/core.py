"""Allocate a process-bigraph core with every Process / Step / emitter
this composite uses registered.

Note that NFSim's wrapper exposes registered link names ``nfsim`` and
``monomer-production`` via its own setup helpers; we register the
classes here under stable names so the document is self-contained.
"""

from process_bigraph import allocate_core
from process_bigraph.emitter import RAMEmitter

from pbg_caspule.processes import CASPULEProcess
from pbg_nfsim.processes import NFSimProcess

from pbg_composite_nfsim_caspule.detector import ObservableDetector


def build_core():
    """Return a freshly allocated core with all classes registered."""
    core = allocate_core()
    core.register_link('CASPULEProcess', CASPULEProcess)
    core.register_link('NFSimProcess', NFSimProcess)
    core.register_link('ObservableDetector', ObservableDetector)
    core.register_link('ram-emitter', RAMEmitter)
    return core
