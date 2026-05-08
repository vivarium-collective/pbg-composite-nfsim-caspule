# pbg-composite-nfsim-caspule

A process-bigraph composite that couples a **spatial bond-aware MD
simulation** ([CASPULE](https://github.com/vivarium-collective/pbg-caspule),
a modified LAMMPS) to a **non-spatial rule-based simulation**
([NFSim/BioNetGen](https://github.com/vivarium-collective/pbg-nfsim))
through a configurable `ObservableDetector` Step.

**[View Interactive Demo Report](https://vivarium-collective.github.io/pbg-composite-nfsim-caspule/)** — four tabbed experiments (decoupled baseline, simple Dimer coupling, NFSim at 5x rate, and a multi-rule polymer-fractionation set-up), each showing the time-series of cross-process coupling, a 3D scatter of CASPULE positions and bond edges, the cluster-size distribution, the detector-rules table, the LAMMPS / BNGL / YAML input files, and a collapsible JSON tree of the full PBG document.

The detector reads CASPULE's live atom positions, atom types, and
global bond list every step, applies a YAML-configured set of
detection rules, and:

1. Removes the matched atoms from CASPULE so the spatial pool no
   longer carries them.
2. Credits each match as one molecule of a named NFSim observable, so
   the non-spatial pool grows at the rate clusters appear in the
   spatial pool.

This is the canonical pattern for handing matter from a *resolved*
representation (atoms + bonds in space) to a *coarse-grained* one
(species counts in a non-spatial pool).

## Architecture

```
                +----------+    positions/bonds/    +----------+    molecule_     +-------+
                |          |    atom_types          |          |    additions     |       |
input.in -----> | CASPULE  | ---------------------> | Detector | ---------------> | NFSim | <-- input.bngl
                |          | <--------------------- |          |                  |       |     (rules)
                +----------+    atoms_to_remove     +----------+                  +-------+
                                                          ^
                                                          | input.yaml (rules)
```

Three input files configure the simulation:

| File | Used by | Purpose |
|---|---|---|
| `*.in` (LAMMPS / CASPULE script) | CASPULE | spatial dynamics, bond formation |
| `*.bngl` (BioNetGen Language) | NFSim | non-spatial reaction rules |
| `*.yaml` (detector rules) | ObservableDetector | which clusters become which molecules |

## Connection table

The full wiring is in [`pbg_composite_nfsim_caspule/wiring.py`](pbg_composite_nfsim_caspule/wiring.py).

| # | Producer | Port | Schema | Case | Consumer | Port |
|---|---|---|---|---|---|---|
| 1 | caspule | positions | `overwrite[list]` | pass-through | detector | positions |
| 2 | caspule | bonds | `overwrite[list]` | pass-through | detector | bonds |
| 3 | caspule | atom_types | `overwrite[list]` | pass-through | detector | atom_types |
| 4 | caspule | num_atoms | `overwrite[integer]` | pass-through | detector | num_atoms |
| 5 | detector | atoms_to_remove | `overwrite[list]` | pass-through | caspule | atoms_to_remove |
| 6 | detector | molecule_additions | `map[float]` (deltas) | pass-through | nfsim | observables (in) |
| 7 | nfsim | observables | `map[float]` (deltas) | self-loop | nfsim | observables (out) |
| 8 | caspule | thermo + cluster scalars | `overwrite[T]` | sink | emitter | — |

Stores 6 and 7 share one path (`stores/species`) so detector additions
and NFSim's own deltas accumulate in a single `map[float]` pool.

## Updates to the wrapped tools

The composite **adds an `atoms_to_remove` input port to
`pbg_caspule.processes.CASPULEProcess`** (`overwrite[list]`). When
non-empty, CASPULE issues a temporary group + `delete_atoms group ...`
before integrating; otherwise the wrapper's behavior is unchanged.
Three new tests in `pbg-caspule/tests/test_processes.py` cover the
removal path (basic deletion, empty no-op, idempotent re-issue).

`pbg-nfsim` is used unchanged.

## Install

```bash
git clone <this repo> pbg-composite-nfsim-caspule
cd pbg-composite-nfsim-caspule

uv venv .venv
source .venv/bin/activate
uv pip install process-bigraph bigraph-schema bigraph-viz pyyaml pytest matplotlib plotly mpich
uv pip install "setuptools<81"          # bionetgen needs pkg_resources
uv pip install -e .

# editable installs of the wrapped tools (sibling clones)
uv pip install -e ../pbg-caspule
uv pip install -e ../pbg-nfsim
```

## Quick start

```python
from process_bigraph import Composite
from pbg_composite_nfsim_caspule import build_core, build_document

doc = build_document(
    caspule_input_file='examples/caspule_dimerize.in',
    nfsim_model_file='examples/nfsim_dimer_kinetics.bngl',
    detector_config='examples/detector_dimer.yaml',
    caspule_interval=0.05,
    nfsim_interval=0.5,
)
sim = Composite({'state': doc}, core=build_core())
sim.run(5.0)
print(sim.state['stores']['species'])    # {'Dimer': ..., 'Active': ...}
print(sim.state['stores']['num_atoms'])  # CASPULE atoms still in the box
```

## Demo

```bash
python demo/demo_report.py
```

Produces `demo/report.html` containing three side-by-side runs:

* **decoupled**: detector outputs routed to dead-end stores. Spatial
  pool stays full, species pool stays at zero — sanity check.
* **coupled**: actual composition. Watch CASPULE atoms drain as
  dimers form, and NFSim's species pool grow.
* **stressed**: same coupling with NFSim ticking 5x faster, so
  Decay drains the pool fast enough to reach quasi-steady state.

## Tests

```bash
pytest
```

Three suites:

* `tests/test_detector.py` — pure-function unit tests on the detector
  rules, no Composite involved.
* `tests/test_assembly.py` — `Composite()` instantiates without
  schema-reconciliation errors, in both coupled and decoupled modes.
* `tests/test_run.py` — end-to-end run; asserts that the species pool
  receives matter *and* that CASPULE loses atoms — the only way both
  can be true is if the full bridge ran.

## Limitations

* The detector currently treats every match as one molecule of one
  named observable. There is no "this cluster produces 2 of A and 1 of B"
  rule; if you need that, layer a second detector or extend
  `DetectorRule`.
* Matched atoms are removed from CASPULE *atomically*; if you need
  some matched atoms to stay (e.g., a chaperone in the cluster
  remains in the spatial pool), add a per-rule `keep_atoms` filter.
* NFSim's BNGL seed-species substitution is regex-based; seed lines
  must use named parameters, not literal numbers, or the substitution
  may bind to the wrong line. See `examples/nfsim_dimer_kinetics.bngl`.
* The composition is one-way at the species layer: NFSim's outputs
  do not flow back into CASPULE as new spatial atoms. Adding that
  reverse coupling would require a CASPULE-side Step that injects
  atoms (LAMMPS `create_atoms` per-step), which is not implemented.
