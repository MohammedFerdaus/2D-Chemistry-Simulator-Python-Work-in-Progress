# Display information
width = 1100
height = 700
sim_width_px = 900
fps = 60
background_color = (15, 15, 20)
simulation_background_color = (20, 20, 30)
panel_background_color = (28, 28, 38)

# Simulation units
simulation_width = 90
simulation_height = 70
scale = 10

# Integration units
dt = 0.01
bond_check_interval = 5
molecule_check_interval = 30

# Thermostat units
gamma = 0.5
temperature = 1.2

# LJ potential units
epsilon = 1.0
lj_cuttoff_factor = 2.5
force_cap = 500

# Element data
elements = {
    'H': {'mass': 1.0, 'max_bonds': 1, 'color': (220, 220, 230), 'sigma': 0.8, 'radius_pixel': 7},
    'O':  {'mass': 16.0, 'max_bonds': 2, 'color': (220, 60, 60), 'sigma': 1.1, 'radius_pixel': 10},
    'N':  {'mass': 14.0, 'max_bonds': 3, 'color': (80, 130, 220), 'sigma': 1.2, 'radius_pixel': 10},
    'He': {'mass': 4.0, 'max_bonds': 0, 'color': (160, 230, 160), 'sigma': 1.0, 'radius_pixel': 8}
}

# Bonding information
bonding_form_factor = 1.4   # r_form = factor * sigma_ij
bonding_break_factor = 2.5  # r_break = factor * r0
spring_k = 300
dissipation = 0.6

bond_lengths = {
    frozenset({'H', 'H'}): 0.74,
    frozenset({'H', 'O'}): 0.96,
    frozenset({'H', 'N'}): 1.01,
    frozenset({'O', 'O'}): 1.21,
    frozenset({'N', 'N'}): 1.10,
    frozenset({'N', 'O'}): 1.15,
}

bond_orders = {            # for double/triple line rendering
    frozenset({'O', 'O'}): 2,
    frozenset({'N', 'N'}): 3,
    frozenset({'N', 'O'}): 2,
}

# Drag information
velocity_scale = 0.05    # drag pixel distance * this = initial sim velocity

# Molecules information
known_molecules = {
    frozenset({('H', 2)}): 'H₂',
    frozenset({('O', 2)}): 'O₂',
    frozenset({('N', 2)}): 'N₂',
    frozenset({('H', 2), ('O', 1)}): 'H₂O',
    frozenset({('H', 3), ('N', 1)}): 'NH₃',
    frozenset({('H', 1), ('O', 1)}): 'OH·',
    frozenset({('H', 2), ('O', 2)}): 'H₂O₂',
    frozenset({('N', 1), ('O', 1)}): 'NO·',
}