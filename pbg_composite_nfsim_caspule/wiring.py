"""Connection table for the CASPULE + NFSim composition.

Single source of truth for cross-process port-to-store wiring. The
document builder reads from this map; tests and the README architecture
diagram should reference it too.

Three classes of stores:

    Spatial (CASPULE -> Detector, pass-through, all overwrite[T]):
        positions, bonds, atom_types, num_atoms

    Action store (Detector -> CASPULE, overwrite[list]):
        atoms_to_remove

    Species pool (Detector -> NFSim and NFSim -> NFSim, bare map[float]
    so deltas compose):
        species

    Sinks (CASPULE outputs that are emitted but not consumed):
        thermo and cluster scalars
"""

WIRING = {
    # spatial pass-through
    ('caspule', 'positions'):           ['stores', 'positions'],
    ('caspule', 'bonds'):               ['stores', 'bonds'],
    ('caspule', 'atom_types'):          ['stores', 'atom_types'],
    ('caspule', 'num_atoms'):           ['stores', 'num_atoms'],
    ('detector', 'positions'):          ['stores', 'positions'],
    ('detector', 'bonds'):               ['stores', 'bonds'],
    ('detector', 'atom_types'):          ['stores', 'atom_types'],
    ('detector', 'num_atoms'):          ['stores', 'num_atoms'],

    # action: detector -> caspule
    ('detector', 'atoms_to_remove'):    ['stores', 'atoms_to_remove'],
    ('caspule',  'atoms_to_remove'):    ['stores', 'atoms_to_remove'],

    # species pool: detector additions + nfsim deltas + nfsim reads
    ('detector', 'molecule_additions'): ['stores', 'species'],
    ('nfsim',    'observables_in'):     ['stores', 'species'],
    ('nfsim',    'observables_out'):    ['stores', 'species'],

    # caspule scalar/thermo sinks
    ('caspule', 'temperature'):         ['stores', 'temperature'],
    ('caspule', 'potential_energy'):    ['stores', 'potential_energy'],
    ('caspule', 'kinetic_energy'):      ['stores', 'kinetic_energy'],
    ('caspule', 'total_energy'):        ['stores', 'total_energy'],
    ('caspule', 'pressure'):            ['stores', 'pressure'],
    ('caspule', 'volume'):              ['stores', 'volume'],
    ('caspule', 'box_dimensions'):      ['stores', 'box_dimensions'],
    ('caspule', 'velocities'):          ['stores', 'velocities'],
    ('caspule', 'num_bonds'):           ['stores', 'num_bonds'],
    ('caspule', 'bonds_by_type'):       ['stores', 'bonds_by_type'],
    ('caspule', 'bond_energy'):         ['stores', 'bond_energy'],
    ('caspule', 'formed_bonds'):        ['stores', 'formed_bonds'],
    ('caspule', 'broken_bonds'):        ['stores', 'broken_bonds'],
    ('caspule', 'num_clusters'):        ['stores', 'num_clusters'],
    ('caspule', 'largest_cluster'):     ['stores', 'largest_cluster'],
    ('caspule', 'cluster_sizes'):       ['stores', 'cluster_sizes'],
}


# Scalar ports the demo emitter records by default.
SCALAR_EMIT_TYPES = {
    'temperature': 'float',
    'potential_energy': 'float',
    'kinetic_energy': 'float',
    'total_energy': 'float',
    'bond_energy': 'float',
    'num_atoms': 'integer',
    'num_bonds': 'integer',
    'formed_bonds': 'integer',
    'broken_bonds': 'integer',
    'num_clusters': 'integer',
    'largest_cluster': 'integer',
}
