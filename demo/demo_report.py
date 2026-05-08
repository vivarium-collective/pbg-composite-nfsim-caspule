"""Demo: pbg-composite-nfsim-caspule cross-process coupling report.

Runs four configurations of the composite and emits a single
self-contained HTML with tabbed navigation. Each tab shows the
PBG document (collapsible), the detector config, the CASPULE input
script, the time-series of cross-process coupling, the final spatial
state (3D scatter + cluster histogram), and a summary table.

Configurations:

    decoupled  -- detector outputs route to dead-end stores. Sanity
                  check: spatial pool stays full, non-spatial empty.
    coupled    -- the simple set-up: a single Dimer rule. Bonded pairs
                  drain into NFSim as Dimer molecules.
    stressed   -- coupled with NFSim ticking 5x faster. Decay drains
                  the pool fast enough to reach quasi-steady state.
    polymer    -- multi-rule detector. Splits the spatial bond graph
                  into Dimer / Trimer / Polymer buckets, only the
                  small ones are removed from CASPULE so the polymer
                  network remains visible in the spatial state view.
"""

import copy
import json
import os

from process_bigraph import Composite, gather_emitter_results

from pbg_composite_nfsim_caspule.core import build_core
from pbg_composite_nfsim_caspule.document import build_document
from pbg_composite_nfsim_caspule.detector import (
    DetectorConfig,
    _parse_detector_config,
    load_detector_config,
)


HERE = os.path.dirname(os.path.abspath(__file__))
EXAMPLES = os.path.normpath(os.path.join(HERE, '..', 'examples'))


# Small but bond-active world used for the dimer-only configurations.
DIMER_CASPULE = open(os.path.join(EXAMPLES, 'caspule_dimerize.in')).read()
DIMER_BNGL    = os.path.join(EXAMPLES, 'nfsim_dimer_kinetics.bngl')
DIMER_RULES   = load_detector_config(
    os.path.join(EXAMPLES, 'detector_dimer.yaml')).rules

# Multi-bond LAMMPS world for the polymer experiment.
POLY_CASPULE  = open(os.path.join(EXAMPLES, 'caspule_polymer.in')).read()
POLY_BNGL     = os.path.join(EXAMPLES, 'nfsim_polymer_kinetics.bngl')
POLY_YAML     = os.path.join(EXAMPLES, 'detector_polymer.yaml')
POLY_RULES    = load_detector_config(POLY_YAML).rules


def _rules_as_dicts(rules):
    """Detector rules as JSON-friendly dicts (the YAML round-trip)."""
    out = []
    for r in rules:
        d = {
            'name': r.name,
            'min_size': r.min_size,
            'remove_from_caspule': r.remove_from_caspule,
        }
        if r.bond_type_cluster is not None:
            d['bond_type_cluster'] = r.bond_type_cluster
        if r.atom_type is not None:
            d['atom_type'] = r.atom_type
        if r.max_size is not None:
            d['max_size'] = r.max_size
        if r.max_bonds is not None:
            d['max_bonds'] = r.max_bonds
        out.append(d)
    return out


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
        'caspule_script': DIMER_CASPULE,
        'nfsim_model':    DIMER_BNGL,
        'detector_rules': _rules_as_dicts(DIMER_RULES),
        'caspule_interval': 0.05,
        'nfsim_interval':   0.5,
        'total_time':       5.0,
        'decoupled':        True,
    },
    {
        'id': 'coupled',
        'title': 'Coupled (simple Dimer rule)',
        'description': (
            'Default wiring with a single rule: bonded pairs of type-1 '
            'atoms (size exactly 2) are pulled out of the spatial pool '
            'and credited as Dimer molecules to NFSim. NFSim then runs '
            'the Dimer <-> Active -> Decay kinetics on the accumulated pool.'
        ),
        'caspule_script': DIMER_CASPULE,
        'nfsim_model':    DIMER_BNGL,
        'detector_rules': _rules_as_dicts(DIMER_RULES),
        'caspule_interval': 0.05,
        'nfsim_interval':   0.5,
        'total_time':       5.0,
        'decoupled':        False,
    },
    {
        'id': 'stressed',
        'title': 'Coupled, NFSim at 5x rate',
        'description': (
            'Same coupling as the simple set-up, but NFSim ticks five '
            'times more often. Decay drains the pool fast enough to '
            'reach quasi-steady state instead of growing monotonically.'
        ),
        'caspule_script': DIMER_CASPULE,
        'nfsim_model':    DIMER_BNGL,
        'detector_rules': _rules_as_dicts(DIMER_RULES),
        'caspule_interval': 0.05,
        'nfsim_interval':   0.1,
        'total_time':       5.0,
        'decoupled':        False,
    },
    {
        'id': 'polymer',
        'title': 'Polymer fractionation (multi-rule)',
        'description': (
            'A more interesting set-up. CASPULE allows up to four bonds '
            'per atom, so atoms grow into branched networks. The '
            'detector applies three rules: Dimer (size 2, removed), '
            'Trimer (size 3, removed) and Polymer (size 4+, observed but '
            'NOT removed). The non-spatial pool fills with size-bucketed '
            'species, while the spatial pool retains the polymer network — '
            'visible in the 3D view at the bottom of this tab.'
        ),
        'caspule_script': POLY_CASPULE,
        'nfsim_model':    POLY_BNGL,
        'detector_rules': _rules_as_dicts(POLY_RULES),
        'caspule_interval': 0.05,
        'nfsim_interval':   0.5,
        'total_time':       5.0,
        'decoupled':        False,
    },
]


# ── helpers ─────────────────────────────────────────────────────────

def _capture_spatial_state(sim):
    s = sim.state.get('stores', {})
    return {
        'num_atoms':     int(s.get('num_atoms') or 0),
        'num_bonds':     int(s.get('num_bonds') or 0),
        'positions':     list(s.get('positions') or []),
        'bonds':         list(s.get('bonds') or []),
        'atom_types':    list(s.get('atom_types') or []),
        'cluster_sizes': list(s.get('cluster_sizes') or []),
        'box_dimensions': list(s.get('box_dimensions') or []),
    }


def _sanitize(obj, depth=0):
    """Recursively convert a doc dict into a JSON-safe structure.

    The Composite document we hand off is already plain dicts/lists/
    strings/numbers, but we deepcopy and round-trip through JSON to be
    sure we never carry process instances or cycles into the report.
    """
    if depth > 12:
        return '... (truncated)'
    if isinstance(obj, dict):
        return {k: _sanitize(v, depth + 1) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_sanitize(v, depth + 1) for v in obj]
    if isinstance(obj, (str, int, float, bool)) or obj is None:
        return obj
    return repr(obj)


def run_config(cfg):
    """Run one config end-to-end. Returns history + doc + spatial snapshot."""
    core = build_core()
    doc = build_document(
        caspule_input_script=cfg['caspule_script'],
        nfsim_model_file=cfg['nfsim_model'],
        detector_rules=cfg['detector_rules'],
        caspule_interval=cfg['caspule_interval'],
        nfsim_interval=cfg['nfsim_interval'],
        nfsim_n_steps=20,
        decoupled=cfg['decoupled'],
        emit_scalars=True,
    )
    doc_for_display = _sanitize(copy.deepcopy(doc))

    sim = Composite({'state': doc}, core=core)

    # Run a small initial slice so the spatial stores have meaningful
    # values before bulk drainage begins. The "initial" snapshot here
    # is post-first-CASPULE-step, which is when the lattice + first
    # bonds become observable.
    warmup = min(cfg['total_time'] * 0.1,
                 cfg['caspule_interval'] * 4)
    sim.run(warmup)
    initial_spatial = _capture_spatial_state(sim)

    sim.run(cfg['total_time'] - warmup)
    final_spatial = _capture_spatial_state(sim)

    history = gather_emitter_results(sim)[('emitter',)]
    return {
        'history': history,
        'doc': doc_for_display,
        'initial_spatial': initial_spatial,
        'final_spatial':  final_spatial,
    }


def history_to_series(history, names=('Dimer', 'Trimer', 'Polymer', 'Active')):
    times = [row.get('time', 0.0) for row in history]
    out = {'time': times}
    for n in names:
        out[n] = [row.get('species', {}).get(n, 0.0) for row in history]
    out['num_atoms'] = [row.get('num_atoms', 0) for row in history]
    out['num_bonds'] = [row.get('num_bonds', 0) for row in history]
    out['formed_bonds'] = [row.get('formed_bonds', 0) for row in history]
    return out


# ── HTML rendering ──────────────────────────────────────────────────

def render_architecture_svg():
    return """
<svg viewBox="0 0 880 360" xmlns="http://www.w3.org/2000/svg" class="arch">
  <defs>
    <marker id="arrow" markerWidth="10" markerHeight="10" refX="9" refY="3"
            orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#5a6"></path>
    </marker>
    <marker id="arrow-r" markerWidth="10" markerHeight="10" refX="9" refY="3"
            orient="auto" markerUnits="strokeWidth">
      <path d="M0,0 L0,6 L9,3 z" fill="#b91c1c"></path>
    </marker>
  </defs>
  <rect x="40"  y="60"  width="200" height="120" rx="10"
        fill="#fef3c7" stroke="#b45309" stroke-width="2"/>
  <text x="140" y="90"  text-anchor="middle" font-weight="700">CASPULE</text>
  <text x="140" y="110" text-anchor="middle" font-size="11" fill="#444">Process</text>
  <text x="140" y="135" text-anchor="middle" font-size="11" fill="#444">bond-aware MD</text>
  <text x="140" y="155" text-anchor="middle" font-size="11" fill="#444">.in input file</text>
  <rect x="340" y="60"  width="200" height="120" rx="10"
        fill="#dbeafe" stroke="#1e40af" stroke-width="2"/>
  <text x="440" y="90"  text-anchor="middle" font-weight="700">ObservableDetector</text>
  <text x="440" y="110" text-anchor="middle" font-size="11" fill="#444">Step</text>
  <text x="440" y="135" text-anchor="middle" font-size="11" fill="#444">cluster -> molecule</text>
  <text x="440" y="155" text-anchor="middle" font-size="11" fill="#444">.yaml input file</text>
  <rect x="640" y="60"  width="200" height="120" rx="10"
        fill="#dcfce7" stroke="#166534" stroke-width="2"/>
  <text x="740" y="90"  text-anchor="middle" font-weight="700">NFSim</text>
  <text x="740" y="110" text-anchor="middle" font-size="11" fill="#444">Process</text>
  <text x="740" y="135" text-anchor="middle" font-size="11" fill="#444">rule-based kinetics</text>
  <text x="740" y="155" text-anchor="middle" font-size="11" fill="#444">.bngl input file</text>
  <line x1="240" y1="100" x2="340" y2="100"
        stroke="#5a6" stroke-width="2" marker-end="url(#arrow)"/>
  <text x="290" y="92" text-anchor="middle" font-size="10">positions, bonds,</text>
  <text x="290" y="105" text-anchor="middle" font-size="10">atom_types</text>
  <line x1="340" y1="160" x2="240" y2="160"
        stroke="#b91c1c" stroke-width="2" marker-end="url(#arrow-r)"/>
  <text x="290" y="178" text-anchor="middle" font-size="10" fill="#7f1d1d">atoms_to_remove</text>
  <line x1="540" y1="120" x2="640" y2="120"
        stroke="#5a6" stroke-width="2" marker-end="url(#arrow)"/>
  <text x="590" y="112" text-anchor="middle" font-size="10">molecule_additions</text>
  <rect x="660" y="240" width="160" height="50" rx="8"
        fill="#fff" stroke="#888" stroke-width="1.5"/>
  <text x="740" y="260" text-anchor="middle" font-size="11" font-weight="700">stores/species</text>
  <text x="740" y="277" text-anchor="middle" font-size="10" fill="#666">map[float], deltas</text>
  <line x1="740" y1="180" x2="740" y2="240"
        stroke="#5a6" stroke-width="1.5" stroke-dasharray="3,3"/>
  <line x1="540" y1="140" x2="660" y2="265"
        stroke="#5a6" stroke-width="1.5" stroke-dasharray="3,3"/>
</svg>
"""


def render_json_tree(obj, key='', open_default=True, depth=0):
    """Render a JSON value as a nested <details> tree.

    Keeps top two levels open by default so the document is immediately
    visible, and lets the user collapse / expand deeper subtrees.
    """
    is_root = depth == 0
    open_attr = ' open' if (open_default and depth < 2) else ''

    def _key_label(k):
        return f'<span class="jt-key">{k}</span>' if k != '' else ''

    if isinstance(obj, dict):
        if not obj:
            return f'{_key_label(key)} <span class="jt-empty">{{}}</span>'
        n = len(obj)
        items = ''.join(
            f'<li>{render_json_tree(v, str(k), open_default, depth + 1)}</li>'
            for k, v in obj.items())
        meta = f'<span class="jt-meta">{{ {n} key{"s" if n != 1 else ""} }}</span>'
        summary = f'{_key_label(key)} {meta}'.strip()
        return (f'<details{open_attr} class="jt jt-obj">'
                f'<summary>{summary}</summary>'
                f'<ul class="jt-list">{items}</ul></details>')

    if isinstance(obj, list):
        n = len(obj)
        if n == 0:
            return f'{_key_label(key)} <span class="jt-empty">[]</span>'
        # Inline short scalar lists.
        if n <= 8 and all(isinstance(v, (int, float, str, bool)) or v is None
                          for v in obj):
            inline = '[' + ', '.join(_render_scalar(v) for v in obj) + ']'
            return f'{_key_label(key)} <span class="jt-inline">{inline}</span>'
        items = ''.join(
            f'<li>{render_json_tree(v, str(i), open_default, depth + 1)}</li>'
            for i, v in enumerate(obj[:30]))
        if n > 30:
            items += f'<li class="jt-more">… {n - 30} more</li>'
        meta = f'<span class="jt-meta">[ {n} item{"s" if n != 1 else ""} ]</span>'
        summary = f'{_key_label(key)} {meta}'.strip()
        return (f'<details{open_attr} class="jt jt-arr">'
                f'<summary>{summary}</summary>'
                f'<ul class="jt-list">{items}</ul></details>')

    return f'{_key_label(key)} {_render_scalar(obj)}'


def _render_scalar(v):
    if v is None:
        return '<span class="jt-null">null</span>'
    if isinstance(v, bool):
        return f'<span class="jt-bool">{str(v).lower()}</span>'
    if isinstance(v, (int, float)):
        return f'<span class="jt-num">{v}</span>'
    if isinstance(v, str):
        # Short strings inline; longer strings get cut for display.
        s = v if len(v) <= 80 else v[:77] + '…'
        return f'<span class="jt-str">"{_html_escape(s)}"</span>'
    return f'<code>{_html_escape(repr(v))}</code>'


def _html_escape(s):
    return (s.replace('&', '&amp;').replace('<', '&lt;')
             .replace('>', '&gt;').replace('"', '&quot;'))


def render_detector_rules(rules):
    if not rules:
        return '<p class="empty">no rules</p>'
    rows = []
    for r in rules:
        match = []
        if 'bond_type_cluster' in r:
            match.append(f'bond_type_cluster=<b>{r["bond_type_cluster"]}</b>')
        if 'atom_type' in r:
            match.append(f'atom_type=<b>{r["atom_type"]}</b>')
        if 'max_bonds' in r:
            match.append(f'max_bonds={r["max_bonds"]}')
        rows.append(
            f'<tr>'
            f'<td><code>{r["name"]}</code></td>'
            f'<td>{" &nbsp; ".join(match)}</td>'
            f'<td>{r["min_size"]}</td>'
            f'<td>{r.get("max_size", "—")}</td>'
            f'<td>{"✓" if r["remove_from_caspule"] else "—"}</td>'
            f'</tr>')
    return ('<table class="rules"><thead><tr>'
            '<th>nfsim observable</th><th>match</th>'
            '<th>min_size</th><th>max_size</th><th>remove?</th>'
            '</tr></thead><tbody>' + ''.join(rows) + '</tbody></table>')


def chart_block(cfg, series, palette):
    """Plotly time-series for the species pool + CASPULE atom count."""
    traces = [
        {'x': series['time'], 'y': series['num_atoms'],
         'type': 'scatter', 'mode': 'lines',
         'name': 'CASPULE atoms', 'yaxis': 'y1',
         'line': {'color': '#b45309', 'width': 2}},
    ]
    for name, color in palette:
        if name in series:
            traces.append({
                'x': series['time'], 'y': series[name],
                'type': 'scatter', 'mode': 'lines',
                'name': f'{name} (NFSim)', 'yaxis': 'y2',
                'line': {'color': color, 'width': 2},
            })
    return {
        'id': f"chart-{cfg['id']}",
        'data': traces,
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


def spatial_3d_block(cfg, snapshot, label):
    """3D scatter of CASPULE atom positions, sampled if very large."""
    pos = snapshot.get('positions') or []
    bonds = snapshot.get('bonds') or []
    sample_cap = 600
    if len(pos) > sample_cap:
        step = max(1, len(pos) // sample_cap)
        pos = pos[::step]
    if not pos:
        return None  # caller will render an empty-state message
    xs = [p[0] for p in pos]
    ys = [p[1] for p in pos]
    zs = [p[2] for p in pos]
    traces = [{
        'x': xs, 'y': ys, 'z': zs,
        'mode': 'markers', 'type': 'scatter3d',
        'name': 'atoms',
        'marker': {'size': 3, 'color': '#b45309', 'opacity': 0.7},
    }]
    # Skip drawing bond edges for large networks — keeps the JSON
    # payload tractable. Up to 250 bond segments is plenty for a
    # human read.
    if 0 < len(bonds) <= 250:
        ex, ey, ez = [], [], []
        for b in bonds:
            _t, a1, a2 = int(b[0]), int(b[1]), int(b[2])
            i, j = a1 - 1, a2 - 1
            if 0 <= i < len(snapshot['positions']) and \
               0 <= j < len(snapshot['positions']):
                p1 = snapshot['positions'][i]
                p2 = snapshot['positions'][j]
                ex += [p1[0], p2[0], None]
                ey += [p1[1], p2[1], None]
                ez += [p1[2], p2[2], None]
        traces.append({
            'x': ex, 'y': ey, 'z': ez,
            'mode': 'lines', 'type': 'scatter3d',
            'name': 'bonds',
            'line': {'color': '#1e40af', 'width': 2},
        })
    return {
        'id': f"scatter-{cfg['id']}-{label}",
        'data': traces,
        'layout': {
            'title': f"spatial state — {label} (n={snapshot['num_atoms']} atoms, "
                     f"{snapshot['num_bonds']} bonds)",
            'margin': {'t': 40, 'r': 0, 'b': 0, 'l': 0},
            'height': 380,
            'scene': {
                'xaxis': {'title': 'x'}, 'yaxis': {'title': 'y'},
                'zaxis': {'title': 'z'},
            },
        },
    }


def cluster_histogram_block(cfg, snapshot, label):
    sizes = snapshot.get('cluster_sizes') or []
    if not sizes:
        return None
    counts = {}
    for s in sizes:
        counts[s] = counts.get(s, 0) + 1
    xs = sorted(counts.keys())
    ys = [counts[s] for s in xs]
    return {
        'id': f"hist-{cfg['id']}-{label}",
        'data': [{
            'x': xs, 'y': ys, 'type': 'bar',
            'marker': {'color': '#1e40af'},
        }],
        'layout': {
            'title': f"cluster-size distribution ({label})",
            'xaxis': {'title': 'cluster size (atoms)'},
            'yaxis': {'title': '# clusters'},
            'margin': {'t': 40, 'r': 30, 'b': 50, 'l': 60},
            'height': 280,
        },
    }


PALETTE_BY_CFG = {
    'decoupled': [('Dimer', '#1e40af'), ('Active', '#16a34a')],
    'coupled':   [('Dimer', '#1e40af'), ('Active', '#16a34a')],
    'stressed':  [('Dimer', '#1e40af'), ('Active', '#16a34a')],
    'polymer':   [('Dimer',   '#1e40af'),
                  ('Trimer',  '#9333ea'),
                  ('Polymer', '#0891b2'),
                  ('Active',  '#16a34a')],
}


def render_panel(cfg, result):
    series = history_to_series(result['history'],
                               names=tuple(n for n, _ in PALETTE_BY_CFG[cfg['id']]))
    main_chart = chart_block(cfg, series, PALETTE_BY_CFG[cfg['id']])
    scatter_initial = spatial_3d_block(cfg, result['initial_spatial'], 'early')
    scatter_final   = spatial_3d_block(cfg, result['final_spatial'],   'end')
    hist_final      = cluster_histogram_block(cfg, result['final_spatial'], 'end')

    extra_charts = [main_chart]
    extra_charts += [c for c in (scatter_initial, scatter_final, hist_final)
                     if c is not None]

    summary = {
        'final time':          series['time'][-1] if series['time'] else None,
        'CASPULE atoms (init)': result['initial_spatial']['num_atoms'],
        'CASPULE atoms (end)':  result['final_spatial']['num_atoms'],
        'CASPULE bonds (end)':  result['final_spatial']['num_bonds'],
        **{f'NFSim {n}': series[n][-1] if series[n] else 0.0
           for n, _ in PALETTE_BY_CFG[cfg['id']]},
    }

    summary_rows = ''.join(
        f'<tr><th>{k}</th><td>{v}</td></tr>'
        for k, v in summary.items()
    )

    detector_rules_html = render_detector_rules(cfg['detector_rules'])
    json_tree_html = render_json_tree(result['doc'], '', open_default=True)

    chart_divs = ''.join(
        f'<div id="{c["id"]}" class="chart"></div>' for c in extra_charts
    )

    return {
        'html': f"""
<section id="tab-{cfg['id']}" class="tab">
  <h2>{cfg['title']}</h2>
  <p>{cfg['description']}</p>

  <div class="meta">
    caspule_interval = {cfg['caspule_interval']} &nbsp;|&nbsp;
    nfsim_interval   = {cfg['nfsim_interval']} &nbsp;|&nbsp;
    total_time       = {cfg['total_time']} &nbsp;|&nbsp;
    decoupled        = {cfg['decoupled']}
  </div>

  <h3>Cross-process time series</h3>
  {chart_divs}

  <h3>Final-run summary</h3>
  <table class="kv">{summary_rows}</table>

  <h3>Observable detector rules</h3>
  {detector_rules_html}

  <details>
    <summary><b>CASPULE input script (.in)</b></summary>
    <pre class="script">{_html_escape(cfg['caspule_script'])}</pre>
  </details>

  <details>
    <summary><b>NFSim model file (.bngl)</b></summary>
    <pre class="script">{_html_escape(open(cfg['nfsim_model']).read())}</pre>
  </details>

  <h3>Process-bigraph document</h3>
  <p class="hint">Click a node to expand or collapse it. Top-level
  process / store nodes are open by default.</p>
  <div class="json-tree">{json_tree_html}</div>

</section>""",
        'charts': extra_charts,
    }


def build_html(panels):
    arch_svg = render_architecture_svg()
    nav = ''.join(
        f'<button data-tab="{p["cfg"]["id"]}">{p["cfg"]["title"]}</button>'
        for p in panels
    )
    sections = ''.join(p['rendered']['html'] for p in panels)

    chart_init_lines = []
    for p in panels:
        for c in p['rendered']['charts']:
            chart_init_lines.append(
                f"Plotly.newPlot('{c['id']}', "
                f"{json.dumps(c['data'])}, {json.dumps(c['layout'])}, "
                "{responsive: true});")
    chart_init = '\n'.join(chart_init_lines)

    css = """
body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
       max-width: 1100px; margin: 2em auto; padding: 0 1em; color: #222; }
h1 { border-bottom: 2px solid #1e40af; padding-bottom: 6px; }
h2 { color: #1e40af; margin-top: 0.4em; }
h3 { color: #444; margin-top: 1.5em; font-size: 14px; text-transform: uppercase;
     letter-spacing: 0.04em; }
.intro { color: #444; line-height: 1.5; }
.arch { width: 100%; max-width: 880px; height: auto; margin: 1em 0;
        border: 1px solid #ddd; border-radius: 6px; padding: 6px; }
nav.tabs { position: sticky; top: 0; z-index: 50; background: #fff;
           border-bottom: 2px solid #1e40af; padding: 8px 0; margin: 1em 0 1.5em; }
nav.tabs button { font: inherit; cursor: pointer; padding: 8px 14px;
                  margin-right: 6px; border: 1px solid #cbd5e1;
                  border-radius: 6px; background: #f1f5f9; color: #1e293b; }
nav.tabs button.active { background: #1e40af; color: #fff; border-color: #1e40af; }
.tab { display: none; }
.tab.active { display: block; }
.meta { font-family: monospace; font-size: 12px; color: #555;
        background: #f8fafc; padding: 8px 12px; border-radius: 4px; }
.chart { width: 100%; margin-bottom: 14px; }
table.rules { border-collapse: collapse; font-size: 13px; margin: 0.5em 0;
              width: 100%; }
table.rules th, table.rules td { border: 1px solid #ddd; padding: 6px 10px;
                                   text-align: left; }
table.rules thead { background: #f1f5f9; }
table.kv { border-collapse: collapse; font-size: 13px; }
table.kv th { text-align: left; padding: 4px 14px 4px 0; color: #555;
              font-weight: 500; }
table.kv td { padding: 4px 0; font-family: monospace; }
pre.script { background: #f1f5f9; padding: 0.7em 1em; border-radius: 6px;
             overflow-x: auto; font-size: 12px; max-height: 320px; }
.hint { color: #666; font-size: 12px; }
.json-tree { font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
             font-size: 12px; line-height: 1.55; background: #f8fafc;
             padding: 12px; border-radius: 6px; max-height: 600px;
             overflow: auto; border: 1px solid #e2e8f0; }
.json-tree details { margin-left: 0; }
.json-tree summary { cursor: pointer; padding: 1px 2px; }
.json-tree summary::marker { color: #94a3b8; }
.json-tree ul.jt-list { list-style: none; padding-left: 1.4em; margin: 0; }
.json-tree li { padding: 1px 0; }
.json-tree .jt-key { color: #1e40af; }
.json-tree .jt-meta { color: #94a3b8; font-size: 11px; }
.json-tree .jt-empty { color: #94a3b8; }
.json-tree .jt-num { color: #b45309; }
.json-tree .jt-bool { color: #b45309; }
.json-tree .jt-str { color: #047857; }
.json-tree .jt-null { color: #94a3b8; font-style: italic; }
.json-tree .jt-inline { color: #047857; }
.json-tree .jt-more { color: #94a3b8; font-style: italic; }
"""
    js = """
function setActive(id) {
  document.querySelectorAll('.tab').forEach(s => s.classList.remove('active'));
  document.querySelectorAll('nav.tabs button').forEach(b => b.classList.remove('active'));
  const sec = document.getElementById('tab-' + id);
  const btn = document.querySelector(`nav.tabs button[data-tab='${id}']`);
  if (sec) sec.classList.add('active');
  if (btn) btn.classList.add('active');
  // Plotly re-layouts charts inside the now-visible section.
  if (sec) {
    sec.querySelectorAll('.chart').forEach(div => {
      if (div.children.length) Plotly.Plots.resize(div);
    });
  }
  if (history && history.replaceState) {
    history.replaceState(null, '', '#' + id);
  }
}
document.addEventListener('DOMContentLoaded', () => {
  document.querySelectorAll('nav.tabs button').forEach(b => {
    b.addEventListener('click', () => setActive(b.dataset.tab));
  });
  const initial = (location.hash || '#' +
    document.querySelector('nav.tabs button').dataset.tab).slice(1);
  setActive(initial);
});
"""

    plotly_cdn = (
        '<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>')

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
through a configurable <code>ObservableDetector</code> Step. Each
simulator is configured by its own input file: a LAMMPS
<code>.in</code> script, a BNGL <code>.bngl</code> model, and a
detector <code>.yaml</code> rule list. Pick an experiment from the
tabs below to inspect its time series, final spatial state, detector
rules, and the full PBG document.
</p>

<h2 style="margin-top: 1.2em;">Architecture</h2>
{arch_svg}

<nav class="tabs">{nav}</nav>

{sections}

<script>
{chart_init}
{js}
</script>

</body></html>"""


def main():
    panels = []
    for cfg in CONFIGS:
        print(f"running config: {cfg['id']} ...")
        result = run_config(cfg)
        rendered = render_panel(cfg, result)
        panels.append({'cfg': cfg, 'result': result, 'rendered': rendered})

    html = build_html(panels)
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
