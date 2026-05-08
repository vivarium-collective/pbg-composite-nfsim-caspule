"""Demo: pbg-composite-nfsim-caspule cross-process coupling report.

Runs three configurations of the composite and emits a single
self-contained HTML showing how the spatial CASPULE simulation hands
matter off to NFSim's non-spatial pool through the ObservableDetector.

Configurations:

    decoupled   -- detector outputs route to dead-end stores. The
                   spatial pool stays full, the non-spatial pool stays
                   empty. Sanity check that nothing flows when the
                   bridge is cut.
    coupled     -- the actual composition. Bonded pairs in CASPULE
                   are pulled out and injected into NFSim as Dimer.
    stressed    -- coupled with NFSim running at a higher rate, so
                   downstream conversion (Dimer -> Active -> Decay)
                   keeps up with the spatial supply.
"""

import json
import os

from process_bigraph import Composite, gather_emitter_results

from pbg_composite_nfsim_caspule.core import build_core
from pbg_composite_nfsim_caspule.document import build_document
from pbg_composite_nfsim_caspule.wiring import WIRING


HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.normpath(os.path.join(HERE, '..', 'examples'))


# A small but bond-active LAMMPS world that runs in seconds.
DEMO_CASPULE_SCRIPT = """
units lj
atom_style bond
boundary p p p
neighbor 0.5 bin
neigh_modify every 1 delay 0 check yes
region box block 0 6 0 6 0 6
create_box 1 box bond/types 1 extra/bond/per/atom 4 extra/special/per/atom 50
mass 1 1.0
lattice sc 1.6
create_atoms 1 box
pair_style soft 2.0
pair_coeff * * 30.0
bond_style harmonic
bond_coeff 1 5.0 1.5
minimize 1e-4 1e-4 100 100
fix nve all nve/limit 0.05
fix lan all langevin 0.4 0.4 1.0 12345
timestep 0.005
fix mkbnd all bond/create 1 1 1 1.8 1 iparam 1 1 jparam 1 1
thermo 5
"""


DEMO_RULES = [
    {'name': 'Dimer', 'bond_type_cluster': 1,
     'min_size': 2, 'max_size': 2, 'remove_from_caspule': True},
]


CONFIGS = [
    {
        'id': 'decoupled',
        'title': 'Decoupled (baseline)',
        'description': (
            'Detector outputs are routed to dead-end stores. Spatial '
            'CASPULE keeps every atom, NFSim never receives matter. '
            'Used as a baseline so the contribution of the bridge is '
            'visible by contrast in the other configurations.'
        ),
        'caspule_interval': 0.05,
        'nfsim_interval': 0.5,
        'total_time': 5.0,
        'decoupled': True,
    },
    {
        'id': 'coupled',
        'title': 'Coupled',
        'description': (
            'Default wiring: ObservableDetector reads CASPULE bonds, '
            'pulls each newly-formed dimer out of the spatial pool '
            '(removing both atoms) and credits it as a Dimer molecule '
            'to NFSim. NFSim then runs Dimer <-> Active -> Decay '
            'kinetics on the accumulated pool.'
        ),
        'caspule_interval': 0.05,
        'nfsim_interval': 0.5,
        'total_time': 5.0,
        'decoupled': False,
    },
    {
        'id': 'stressed',
        'title': 'Coupled — NFSim at 5x rate',
        'description': (
            'Same coupling, but NFSim ticks five times more often. '
            'Active forms accumulate faster and Decay drains the pool '
            'in real time, so the species pool reaches a quasi-steady '
            'state instead of growing monotonically.'
        ),
        'caspule_interval': 0.05,
        'nfsim_interval': 0.1,
        'total_time': 5.0,
        'decoupled': False,
    },
]


def run_config(cfg):
    """Run one config end-to-end and return (history, doc).

    `history` is the RAM emitter's results: a list of timestamped
    snapshots of the scalar stores plus species pool. `doc` is the
    Composite document used, returned for display only.
    """
    core = build_core()
    doc = build_document(
        caspule_input_script=DEMO_CASPULE_SCRIPT,
        nfsim_model_file=os.path.join(EXAMPLES, 'nfsim_dimer_kinetics.bngl'),
        detector_rules=DEMO_RULES,
        caspule_interval=cfg['caspule_interval'],
        nfsim_interval=cfg['nfsim_interval'],
        nfsim_n_steps=20,
        decoupled=cfg['decoupled'],
        emit_scalars=True,
    )
    sim = Composite({'state': doc}, core=core)
    sim.run(cfg['total_time'])
    history = gather_emitter_results(sim)[('emitter',)]
    return history, doc


def history_to_series(history):
    """Pull a small set of named series out of the emitter dump."""
    times = [row.get('time', 0.0) for row in history]
    num_atoms = [row.get('num_atoms', 0) for row in history]
    num_bonds = [row.get('num_bonds', 0) for row in history]
    formed = [row.get('formed_bonds', 0) for row in history]
    dimer = [row.get('species', {}).get('Dimer', 0.0) for row in history]
    active = [row.get('species', {}).get('Active', 0.0) for row in history]
    return {
        'time': times,
        'num_atoms': num_atoms,
        'num_bonds': num_bonds,
        'formed_bonds': formed,
        'dimer': dimer,
        'active': active,
    }


def render_architecture_svg():
    """Hand-laid box-and-arrow SVG of the composition."""
    return """
<svg viewBox="0 0 880 360" xmlns="http://www.w3.org/2000/svg" class="arch">
  <defs>
    <marker id="arrow" markerWidth="10" markerHeight="10" refX="9" refY="3"
            orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#5a6"></path>
    </marker>
  </defs>

  <!-- CASPULE -->
  <rect x="40"  y="60"  width="200" height="120" rx="10"
        fill="#fef3c7" stroke="#b45309" stroke-width="2"/>
  <text x="140" y="90"  text-anchor="middle" font-weight="700">CASPULE</text>
  <text x="140" y="110" text-anchor="middle" font-size="11" fill="#444">Process</text>
  <text x="140" y="135" text-anchor="middle" font-size="11" fill="#444">bond-aware MD</text>
  <text x="140" y="155" text-anchor="middle" font-size="11" fill="#444">.in input file</text>

  <!-- Detector -->
  <rect x="340" y="60"  width="200" height="120" rx="10"
        fill="#dbeafe" stroke="#1e40af" stroke-width="2"/>
  <text x="440" y="90"  text-anchor="middle" font-weight="700">ObservableDetector</text>
  <text x="440" y="110" text-anchor="middle" font-size="11" fill="#444">Step</text>
  <text x="440" y="135" text-anchor="middle" font-size="11" fill="#444">cluster -> molecule</text>
  <text x="440" y="155" text-anchor="middle" font-size="11" fill="#444">.yaml input file</text>

  <!-- NFSim -->
  <rect x="640" y="60"  width="200" height="120" rx="10"
        fill="#dcfce7" stroke="#166534" stroke-width="2"/>
  <text x="740" y="90"  text-anchor="middle" font-weight="700">NFSim</text>
  <text x="740" y="110" text-anchor="middle" font-size="11" fill="#444">Process</text>
  <text x="740" y="135" text-anchor="middle" font-size="11" fill="#444">rule-based kinetics</text>
  <text x="740" y="155" text-anchor="middle" font-size="11" fill="#444">.bngl input file</text>

  <!-- positions/bonds/atom_types -->
  <line x1="240" y1="100" x2="340" y2="100"
        stroke="#5a6" stroke-width="2" marker-end="url(#arrow)"/>
  <text x="290" y="92" text-anchor="middle" font-size="10">positions, bonds,</text>
  <text x="290" y="105" text-anchor="middle" font-size="10">atom_types</text>

  <!-- atoms_to_remove (back to CASPULE) -->
  <line x1="340" y1="160" x2="240" y2="160"
        stroke="#b91c1c" stroke-width="2" marker-end="url(#arrow)"/>
  <text x="290" y="178" text-anchor="middle" font-size="10" fill="#7f1d1d">atoms_to_remove</text>

  <!-- molecule_additions -->
  <line x1="540" y1="120" x2="640" y2="120"
        stroke="#5a6" stroke-width="2" marker-end="url(#arrow)"/>
  <text x="590" y="112" text-anchor="middle" font-size="10">molecule_additions</text>

  <!-- species (NFSim self) -->
  <path d="M 740 200 Q 740 250 740 200" stroke="#16a34a" stroke-width="2"
        fill="none" marker-end="url(#arrow)"/>
  <text x="740" y="230" text-anchor="middle" font-size="10">species (delta)</text>

  <!-- shared store -->
  <rect x="660" y="240" width="160" height="50" rx="8"
        fill="#fff" stroke="#888" stroke-width="1.5"/>
  <text x="740" y="260" text-anchor="middle" font-size="11" font-weight="700">stores/species</text>
  <text x="740" y="277" text-anchor="middle" font-size="10" fill="#666">map[float], deltas</text>

  <line x1="740" y1="180" x2="740" y2="240"
        stroke="#5a6" stroke-width="1.5" stroke-dasharray="3,3"/>
  <line x1="540" y1="140" x2="660" y2="265"
        stroke="#5a6" stroke-width="1.5" stroke-dasharray="3,3"/>

  <text x="40"  y="330" font-size="10" fill="#666">
    spatial pool (LAMMPS atoms / bonds)
  </text>
  <text x="640" y="330" font-size="10" fill="#666">
    non-spatial species pool
  </text>
</svg>
"""


def chart_block(cfg, series):
    """Inline Plotly chart for a single config."""
    return {
        'id': f"chart-{cfg['id']}",
        'data': [
            {'x': series['time'], 'y': series['num_atoms'],
             'type': 'scatter', 'mode': 'lines',
             'name': 'CASPULE atoms',
             'yaxis': 'y1',
             'line': {'color': '#b45309', 'width': 2}},
            {'x': series['time'], 'y': series['dimer'],
             'type': 'scatter', 'mode': 'lines',
             'name': 'Dimer (NFSim)',
             'yaxis': 'y2',
             'line': {'color': '#1e40af', 'width': 2}},
            {'x': series['time'], 'y': series['active'],
             'type': 'scatter', 'mode': 'lines',
             'name': 'Active (NFSim)',
             'yaxis': 'y2',
             'line': {'color': '#16a34a', 'width': 2}},
        ],
        'layout': {
            'title': cfg['title'],
            'xaxis': {'title': 'time (lj units)'},
            'yaxis':  {'title': 'CASPULE atoms', 'side': 'left',
                       'titlefont': {'color': '#b45309'},
                       'tickfont':  {'color': '#b45309'}},
            'yaxis2': {'title': 'NFSim count', 'side': 'right',
                       'overlaying': 'y',
                       'titlefont': {'color': '#1e40af'},
                       'tickfont':  {'color': '#1e40af'}},
            'legend': {'orientation': 'h'},
            'margin': {'t': 50, 'r': 60, 'b': 50, 'l': 60},
            'height': 360,
        },
    }


def build_html(panels, arch_svg):
    """Single-file HTML report with embedded Plotly charts."""
    plotly_cdn = (
        '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>')
    chart_specs = [chart_block(p['cfg'], p['series']) for p in panels]
    chart_init = '\n'.join(
        f"Plotly.newPlot('{c['id']}', "
        f"{json.dumps(c['data'])}, {json.dumps(c['layout'])}, "
        "{responsive: true});"
        for c in chart_specs
    )

    panel_html = []
    for p, c in zip(panels, chart_specs):
        cfg = p['cfg']
        panel_html.append(f"""
<section class="panel">
  <h2>{cfg['title']}</h2>
  <p>{cfg['description']}</p>
  <div class="meta">
    caspule_interval = {cfg['caspule_interval']} &nbsp;|&nbsp;
    nfsim_interval   = {cfg['nfsim_interval']} &nbsp;|&nbsp;
    total_time       = {cfg['total_time']} &nbsp;|&nbsp;
    decoupled        = {cfg['decoupled']}
  </div>
  <div id="{c['id']}" class="chart"></div>
  <details>
    <summary>final state</summary>
    <pre>{json.dumps(p['final'], indent=2)}</pre>
  </details>
</section>
""")

    css = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }
h1 { border-bottom: 2px solid #1e40af; padding-bottom: 6px; }
h2 { color: #1e40af; margin-top: 2em; }
.arch { width: 100%; max-width: 880px; height: auto; margin: 1em 0; }
.panel { border: 1px solid #ddd; padding: 1em 1.5em; border-radius: 8px;
         margin-bottom: 2em; background: #fafafa; }
.meta { font-family: monospace; font-size: 12px; color: #555; }
.chart { width: 100%; }
pre { background: #f1f5f9; padding: 0.5em 1em; border-radius: 6px;
      overflow-x: auto; font-size: 12px; }
.intro { color: #444; line-height: 1.5; }
"""

    return f"""<!doctype html>
<html><head>
<meta charset="utf-8">
<title>pbg-composite-nfsim-caspule report</title>
{plotly_cdn}
<style>{css}</style>
</head><body>

<h1>pbg-composite-nfsim-caspule</h1>
<p class="intro">
A process-bigraph composite that hands matter from a spatial bond-aware
LAMMPS / CASPULE simulation into NFSim's non-spatial rule-based pool
through a configurable <code>ObservableDetector</code> Step. The
detector reads the live bond network, treats user-specified spatial
patterns as the formation event for a named NFSim molecule, and
removes the matched atoms from CASPULE so they are not double-counted.
Each simulator is configured by its own input file: a LAMMPS
<code>.in</code> script, a BNGL <code>.bngl</code> model, and a
detector <code>.yaml</code> rule list.
</p>

<h2>Architecture</h2>
{arch_svg}

{''.join(panel_html)}

<script>
{chart_init}
</script>

</body></html>"""


def main():
    panels = []
    for cfg in CONFIGS:
        print(f"running config: {cfg['id']} ...")
        history, _doc = run_config(cfg)
        series = history_to_series(history)
        final = {
            'time': series['time'][-1] if series['time'] else None,
            'num_atoms': series['num_atoms'][-1] if series['num_atoms'] else None,
            'num_bonds': series['num_bonds'][-1] if series['num_bonds'] else None,
            'Dimer': series['dimer'][-1] if series['dimer'] else None,
            'Active': series['active'][-1] if series['active'] else None,
        }
        panels.append({'cfg': cfg, 'series': series, 'final': final})

    arch_svg = render_architecture_svg()
    html = build_html(panels, arch_svg)

    output_path = os.path.join(HERE, 'report.html')
    with open(output_path, 'w') as f:
        f.write(html)
    print(f"wrote {output_path}")

    try:
        import webbrowser
        webbrowser.open('file://' + output_path)
    except Exception:
        pass


if __name__ == '__main__':
    main()
