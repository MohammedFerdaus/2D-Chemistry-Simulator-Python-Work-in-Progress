from collections import Counter, deque
from config import *
from physics import *

class BondManager:
    def __init__(self):
        self.bonds = set()
        self.flash_events = []

    def update_bonds(self, atoms):
        # Filter out held atoms
        free_atoms = [atom for atom in atoms if not atom.held]

        for i in range(len(free_atoms)):
            for j in range(i + 1, len(free_atoms)):
                atom_i = free_atoms[i]
                atom_j = free_atoms[j]

                # Skip if either atom has no open valence slots
                if (len(atom_i.bonds) >= atom_i.max_bonds or
                    len(atom_j.bonds) >= atom_j.max_bonds):
                    continue

                # Skip if pair is already bonded
                if frozenset({atom_i.id, atom_j.id}) in self.bonds:
                    continue

                # Compute bond formation threshold
                sigma_ij = (atom_i.sigma + atom_j.sigma) / 2
                r_form = sigma_ij * config.bonding_form_factor

                # Compute distance
                r = np.linalg.norm(atom_j.pos - atom_i.pos)

                # Form bond if close enough
                if r < r_form:
                    self.form_bond(atom_i, atom_j)

    def form_bond(self, atom_i, atom_j):
        # Mutual bond registration
        atom_i.bonds.append(atom_j)
        atom_j.bonds.append(atom_i)
        self.bonds.add(frozenset({atom_i.id, atom_j.id}))

        # Mass weighted energy dissipation
        v_rel = atom_i.vel - atom_j.vel
        total_mass = atom_i.mass + atom_j.mass
        delta_v = config.dissipation * v_rel

        atom_i.vel -= (atom_j.mass / total_mass) * delta_v
        atom_j.vel += (atom_i.mass / total_mass) * delta_v

        # Flash event at bond midpoint
        midpoint = (atom_i.pos + atom_j.pos) / 2
        self.flash_events.append((midpoint.copy(), 0))

    def break_bond(self, atom_i, atom_j):
        # Remove atom 
        atom_i.bonds.remove(atom_j)
        atom_j.bonds.remove(atom_i)
        self.bonds.discard(frozenset({atom_i.id, atom_j.id}))

def detect_molecules(atoms):
    visited = set()
    molecule_counts = Counter()

    for atom in atoms:
        if atom.id in visited:
            continue

        # BFS to collect connected component
        component = []
        queue = deque([atom])

        while queue:
            current = queue.popleft()

            if current.id in visited:
                continue

            visited.add(current.id)
            component.append(current)

            for neighbour in current.bonds:
                if neighbour.id not in visited:
                    queue.append(neighbour)

        # Count elements in this molecule
        counts = Counter(a.element for a in component)

        # Create canonical key
        key = frozenset(counts.items())

        # Convert to readable chemical formula
        name = known_molecules.get(key, 'unknown')

        molecule_counts[name] += 1

    return molecule_counts
