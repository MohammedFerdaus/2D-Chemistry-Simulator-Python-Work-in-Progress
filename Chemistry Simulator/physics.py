import numpy as np
import itertools
import config

id_counter = itertools.count()

class Atom:
    def __init__(self, element, position):
        # Identity
        self.id = next(id_counter)
        self.element = element

        # Kinematics
        self.pos = position.copy()
        self.vel = np.zeros(2, dtype=np.float64)
        self.force = np.zeros(2, dtype=np.float64)
        self.force_old = np.zeros(2, dtype=np.float64)

        # Properties from config
        self.mass = config.elements[element]['mass']
        self.sigma = config.elements[element]['sigma']
        self.max_bonds = config.elements[element]['max_bonds']

        # State
        self.bonds = []
        self.held = True


def lj_force(atom_i, atom_j):
    # Displacement vector and distance
    d = atom_j.pos - atom_i.pos
    r = max(np.linalg.norm(d), 0.1)

    # Mixed sigma and cutoff check
    sigma_ij = (atom_i.sigma + atom_j.sigma) / 2
    if r >= config.lj_cutoff_factor * sigma_ij:
        return

    # Unit vector
    u = d / r

    # LJ force scalar
    s = sigma_ij / r
    F = (24 * config.epsilon / r) * (2 * s**12 - s**6)
    F = np.clip(F, -config.force_cap, config.force_cap)

    # Apply Newton's third law
    F_vec = F * u
    atom_i.force += F_vec
    atom_j.force -= F_vec


def spring_force(atom_i, atom_j, bond_manager):
    # Displacement vector and distance
    d = atom_j.pos - atom_i.pos
    r = max(np.linalg.norm(d), 0.1)
    u = d / r

    # Equilibrium bond length from config
    key = frozenset({atom_i.element, atom_j.element})
    r0 = config.bond_lengths[key]

    # Break bond if stretched past threshold
    if r >= config.bonding_break_factor * r0:
        bond_manager.break_bond(atom_i, atom_j)
        return

    # Harmonic spring force
    F_vec = config.spring_k * (r - r0) * u
    atom_i.force += F_vec
    atom_j.force -= F_vec


def langevin_force(atom):
    # Drag term — opposes velocity
    F_drag = -config.gamma * atom.mass * atom.vel

    # Random kick — fluctuation-dissipation theorem prefactor
    sigma_noise = np.sqrt(2 * config.gamma * atom.mass * config.temperature / config.dt)
    F_rand = sigma_noise * np.random.normal(0, 1, size=2)

    atom.force += F_drag + F_rand


def angle_force(atom):
    # Only applies to atoms with 2 or more bonds
    if len(atom.bonds) < 2:
        return

    # Only applies to elements with a defined target angle
    if atom.element not in config.bond_angles:
        return

    theta_0 = np.radians(config.bond_angles[atom.element])

    # Loop over unique pairs of bonded neighbours
    for i in range(len(atom.bonds)):
        for j in range(i + 1, len(atom.bonds)):
            n1 = atom.bonds[i]
            n2 = atom.bonds[j]

            # Vectors from centre atom to each neighbour
            r1 = n1.pos - atom.pos
            r2 = n2.pos - atom.pos

            d1 = max(np.linalg.norm(r1), 0.1)
            d2 = max(np.linalg.norm(r2), 0.1)

            r1_hat = r1 / d1
            r2_hat = r2 / d2

            # Current angle — clamped to avoid acos domain errors
            cos_theta = np.clip(np.dot(r1_hat, r2_hat), -1.0, 1.0)
            theta = np.arccos(cos_theta)

            # Restoring torque toward target angle
            tau = -config.angle_spring_k * (theta - theta_0)

            # Perpendicular force directions
            f1_dir = (r2_hat - cos_theta * r1_hat) / d1
            f2_dir = (r1_hat - cos_theta * r2_hat) / d2

            f1 = tau * f1_dir
            f2 = tau * f2_dir

            # Apply to neighbours, reaction force to centre atom
            n1.force += f1
            n2.force += f2
            atom.force -= (f1 + f2)


def apply_boundary(atom):
    # Atom radius in sim units
    radius_sim = config.elements[atom.element]['radius_pixel'] / config.scale

    # x-axis boundaries
    lo_x = radius_sim
    hi_x = config.simulation_width - radius_sim

    if atom.pos[0] < lo_x:
        atom.pos[0] = lo_x
        atom.vel[0] = abs(atom.vel[0])

    if atom.pos[0] > hi_x:
        atom.pos[0] = hi_x
        atom.vel[0] = -abs(atom.vel[0])

    # y-axis boundaries
    lo_y = radius_sim
    hi_y = config.simulation_height - radius_sim

    if atom.pos[1] < lo_y:
        atom.pos[1] = lo_y
        atom.vel[1] = abs(atom.vel[1])

    if atom.pos[1] > hi_y:
        atom.pos[1] = hi_y
        atom.vel[1] = -abs(atom.vel[1])


def velocity_verlet(atoms, bond_manager, dt):
    # Filter out held atoms
    free_atoms = [atom for atom in atoms if not atom.held]

    # Position update using old acceleration
    for atom in free_atoms:
        a_old = atom.force_old / atom.mass
        atom.pos += atom.vel * dt + 0.5 * a_old * dt**2

    # Boundary conditions after position update
    for atom in free_atoms:
        apply_boundary(atom)

    # Store old forces and zero accumulators
    for atom in free_atoms:
        atom.force_old = atom.force.copy()
        atom.force = np.zeros(2, dtype=np.float64)

    # Pairwise forces — bonded get spring, unbound get LJ
    for i in range(len(free_atoms)):
        for j in range(i + 1, len(free_atoms)):
            atom_i = free_atoms[i]
            atom_j = free_atoms[j]

            if frozenset({atom_i.id, atom_j.id}) in bond_manager.bonds:
                spring_force(atom_i, atom_j, bond_manager)
            else:
                lj_force(atom_i, atom_j)

    # Angle forces — enforces molecular geometry
    for atom in free_atoms:
        angle_force(atom)

    # Langevin thermostat — maintains temperature
    for atom in free_atoms:
        langevin_force(atom)

    # Velocity update using average of old and new acceleration
    for atom in free_atoms:
        a_old = atom.force_old / atom.mass
        a_new = atom.force / atom.mass
        atom.vel += 0.5 * (a_old + a_new) * dt
