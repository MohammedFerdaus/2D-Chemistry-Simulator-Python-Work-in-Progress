import numpy as np
import itertools
import config

id_counter = itertools.count()

# Atom class
class Atom:
    def __init__(self, element, position):
        self.id = next(id_counter)
        self.element = element
        self.pos = position.copy()

        self.ve = np.zeros(2, dtype=np.float64)
        self.force = np.zeros(2, dtype=np.float64)
        self.force_old = np.zeros(2, dtype=np.float64)

        self.mass = config.elements[element]['mass']
        self.sigma = config.elements[element]['sigma']
        self.max_bonds = config.elements[element]['max_bonds']
        
        self.bonds = []
        self.held = True

# Lennard-Jones force function
def lj_force(atom_i, atom_j):
    # Displacement vector
    d = atom_j.pos - atom_i.pos

    # Scalar distance with overlap guard
    r = np.linalg.norm(d)
    r = max(r, 0.1)

    # Cutoff
    sigma_ij = (atom_i.sigma + atom_j.sigma) / 2
    cutoff = config.lj_cutoff_factor * sigma_ij
    if r >= cutoff:
        return

    # Unit vector
    u = d / r

    # Force scalar
    s = sigma_ij / r
    F = (24 * config.epsilon / r) * (2 * s**12 - s**6)

    # Force cap
    F = np.clip(F, -config.force_cap, config.force_cap)    
    
    # Apply to both atoms
    F_vec = F * u
    atom_i.force += F_vec
    atom_j.force -= F_vec

# Spring force function
def spring_force(atom_i, atom_j, bond_manager):
    # Displacement and distance
    d = atom_j.pos - atom_i.pos
    r = np.linalg.norm(d)
    r = max(r, 0.1)
    u = d/r

    # Equilibrium length
    key = frozenset({atom_i.element, atom_j.element})
    r0 = config.bond_lengths[key]

    # Bond breaking
    r_break = config.bonding_break_factor * r0
    if r >= r_break:
        bond_manager.break_bond(atom_i, atom_j)
        return
    
    # Spring scalar force
    F = config.spring_k * (r - r0)

    # Apply to both atoms
    F_vec = F * u
    atom_i.force += F_vec
    atom_j.force -= F_vec

# Langevin force function
def langevin_force(atom):
    # Drag term
    F_drag = (-1 * config.gamma) * atom.mass * atom.vel

    # Noise prefactor
    sigma_noise = np.sqrt(2 * config.gamma * atom.mass * config.temperature / config.dt)

    # Random kick
    eta = np.random.normal(0, 1, size=2)
    F_rand = sigma_noise * eta

    # Apply both
    atom.force += F_drag + F_rand

# Boundary function
def apply_boundary(atom):
    # get atom radius in sim units
    radius_px = config.elements[atom.element]['radius_pixel']
    radius_sim = radius_px / config.scale

    # x-axis
    lo_x = radius_sim
    hi_x = config.simulation_width - radius_sim

    if atom.pos[0] < lo_x:
        atom.pos[0] = lo_x
        atom.vel[0] = abs(atom.vel[0])

    if atom.pos[0] > hi_x:
        atom.pos[0] = hi_x
        atom.vel[0] = -abs(atom.vel[0])

    # y-axis
    lo_y = radius_sim
    hi_y = config.simulation_height - radius_sim

    if atom.pos[1] < lo_y:
        atom.pos[1] = lo_y
        atom.vel[1] = abs(atom.vel[1])

    if atom.pos[1] > hi_y:
        atom.pos[1] = hi_y
        atom.vel[1] = abs(atom.vel[1])

# Velocity verlet function
def velocity_verlet(atoms, bond_manager, dt):

    # Filter out held atoms once and reuse the list
    free_atoms = [atom for atom in atoms if not atom.held]

    # Position update, store old acceleration
    for atom in free_atoms:
        a_old = atom.force_old / atom.mass
        atom.pos += atom.vel * dt + 0.5 * a_old * dt**2

    # Boundary conditions
    for atom in free_atoms:
        apply_boundary(atom)

    # Zero force accumulators and stold old atom force
    for atom in free_atoms:
        atom.force_old = atom.force.copy()
        atom.force = np.zeros(2, dtype=np.float64)

    # Pairwise forces
    for i in range(len(free_atoms)):
        for j in range(i + 1, len(free_atoms)):
            atom_i = free_atoms[i]
            atom_j = free_atoms[j]

            bonded = frozenset({atom_i.id, atom_j.id} in bond_manager.bonds)

            if bonded:
                spring_force(atom_i, atom_j, bond_manager)
            else:
                lj_force(atom_i, atom_j)

    # Langevin per atom
    for atom in free_atoms:
        langevin_force(atom)

    # Velocity update
    for atom in free_atoms:
        a_old = atom.force_old / atom.mass
        a_new = atom.force / atom.mass
        atom.vel += 0.5 * (a_old + a_new) * dt