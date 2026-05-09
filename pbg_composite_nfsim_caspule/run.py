"""Run a pbg-composite-nfsim-caspule simulation from a JSON document.

The PBG JSON file is a Composite document — a Python dict serialised
as JSON — that names the three Process / Step instances (CASPULE,
ObservableDetector, NFSim), their wiring, the shared stores, and the
emitter. Each process's input-file path inside the JSON
(``caspule.config.input_file``, ``detector.config.config_file``,
``nfsim.config.model_file``) is resolved **relative to the JSON
file's directory**, so the bundled examples work without any path
gymnastics.

CLI:

    python -m pbg_composite_nfsim_caspule examples/composite_dimer.pbg.json
    python -m pbg_composite_nfsim_caspule examples/composite_dimer.pbg.json --total-time 2 --output run.json

Override any of the three input files at run time:

    python -m pbg_composite_nfsim_caspule composite.pbg.json \\
        --caspule my.in --detector my.yaml --nfsim my.bngl
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import Any

from process_bigraph import Composite, gather_emitter_results

from pbg_composite_nfsim_caspule.core import build_core


# Process key in the document → config field that holds an input-file path.
INPUT_FILE_FIELDS = {
    'caspule':  ('input_file', 'input_script'),
    'detector': ('config_file',),
    'nfsim':    ('model_file',),
}


def load_document(pbg_path: str) -> tuple[dict, str]:
    """Load a PBG JSON file. Returns (document, base_dir).

    ``base_dir`` is the directory the JSON lives in — used for resolving
    relative input-file paths in the same step.
    """
    pbg_path = os.path.abspath(pbg_path)
    base_dir = os.path.dirname(pbg_path)
    with open(pbg_path) as f:
        doc = json.load(f)
    return doc, base_dir


def resolve_input_paths(doc: dict, base_dir: str) -> None:
    """Mutate ``doc`` so all input-file paths are absolute.

    Relative paths are resolved against ``base_dir``. Empty strings and
    paths that already point at an existing absolute path are left
    alone. Logical fields like ``input_script`` (inline LAMMPS source)
    are also passed through unchanged — they aren't paths.
    """
    for proc_key in ('caspule', 'detector', 'nfsim'):
        if proc_key not in doc:
            continue
        cfg = doc[proc_key].get('config') or {}
        for field in INPUT_FILE_FIELDS.get(proc_key, ()):
            if field == 'input_script':
                continue  # inline source, not a path
            v = cfg.get(field)
            if not v or os.path.isabs(v):
                continue
            cfg[field] = os.path.abspath(os.path.join(base_dir, v))


def apply_overrides(doc: dict, overrides: dict[str, str]) -> None:
    """Apply CLI --caspule / --detector / --nfsim path overrides."""
    map_ = {
        'caspule':  ('caspule',  'input_file'),
        'detector': ('detector', 'config_file'),
        'nfsim':    ('nfsim',    'model_file'),
    }
    for cli_key, path in overrides.items():
        if not path:
            continue
        proc_key, field = map_[cli_key]
        if proc_key not in doc:
            raise ValueError(
                f'cannot override {cli_key}: no "{proc_key}" entry in document')
        doc[proc_key].setdefault('config', {})
        doc[proc_key]['config'][field] = os.path.abspath(path)
        # If the user gave an explicit input_file, clear the inline
        # input_script so it doesn't shadow.
        if cli_key == 'caspule':
            doc['caspule']['config']['input_script'] = ''


def _final_state_summary(sim) -> dict[str, Any]:
    s = sim.state.get('stores', {})
    return {
        'global_time':       sim.state.get('global_time'),
        'caspule_atoms':     int(s.get('num_atoms') or 0),
        'caspule_bonds':     int(s.get('num_bonds') or 0),
        'caspule_clusters':  int(s.get('num_clusters') or 0),
        'caspule_largest':   int(s.get('largest_cluster') or 0),
        'nfsim_species':     dict(s.get('species') or {}),
    }


def run(pbg_path: str, *,
        caspule: str | None = None,
        detector: str | None = None,
        nfsim: str | None = None,
        total_time: float = 5.0,
        output: str | None = None,
        quiet: bool = False) -> dict[str, Any]:
    """Programmatic entry-point. Same semantics as the CLI."""
    doc, base_dir = load_document(pbg_path)
    resolve_input_paths(doc, base_dir)
    apply_overrides(doc, {'caspule': caspule, 'detector': detector,
                          'nfsim': nfsim})

    core = build_core()
    sim = Composite({'state': doc}, core=core)
    sim.run(total_time)

    summary = _final_state_summary(sim)

    if not quiet:
        print(f"final time:           {summary['global_time']}")
        print(f"CASPULE atoms (end):  {summary['caspule_atoms']}")
        print(f"CASPULE bonds (end):  {summary['caspule_bonds']}")
        print(f"CASPULE clusters:     {summary['caspule_clusters']}"
              f" (largest {summary['caspule_largest']})")
        print(f"NFSim species:        {summary['nfsim_species']}")

    if output:
        history = gather_emitter_results(sim).get(('emitter',), [])
        with open(output, 'w') as f:
            json.dump({'summary': summary, 'history': history},
                      f, indent=2, default=str)
        if not quiet:
            print(f"wrote {len(history)} emitter snapshots to {output}")

    return summary


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog='pbg-composite-nfsim-caspule',
        description=(
            'Run a CASPULE + ObservableDetector + NFSim composite from a '
            'PBG JSON document. Bundled examples live under examples/ in '
            'this repo.'),
    )
    parser.add_argument(
        'pbg_file',
        help='Path to the PBG document JSON (e.g. examples/composite_dimer.pbg.json)')
    parser.add_argument(
        '--caspule', metavar='PATH',
        help='Override CASPULE input file (.in / LAMMPS script)')
    parser.add_argument(
        '--detector', metavar='PATH',
        help='Override ObservableDetector config (.yaml)')
    parser.add_argument(
        '--nfsim', metavar='PATH',
        help='Override NFSim model file (.bngl)')
    parser.add_argument(
        '--total-time', type=float, default=5.0,
        help='Total simulation time in CASPULE\'s lj-time units (default 5.0)')
    parser.add_argument(
        '--output', metavar='PATH',
        help='Write final summary + emitter history as JSON to this file')
    parser.add_argument(
        '--quiet', action='store_true',
        help='Suppress stdout summary')
    args = parser.parse_args(argv)

    run(
        args.pbg_file,
        caspule=args.caspule,
        detector=args.detector,
        nfsim=args.nfsim,
        total_time=args.total_time,
        output=args.output,
        quiet=args.quiet,
    )
    return 0


if __name__ == '__main__':
    sys.exit(main())
