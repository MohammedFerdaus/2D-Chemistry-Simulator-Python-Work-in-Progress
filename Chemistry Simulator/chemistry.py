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

                # Bond formation threshold using Lorentz-Berthelot mixing
                sigma_ij = (atom_i.sigma + atom_j.sigma) / 2
                r_form = sigma_ij * bonding_form_factor

                # Compute distance and form bond if close enough
                r = np.linalg.norm(atom_j.pos - atom_i.pos)
                if r < r_form:
                    self.form_bond(atom_i, atom_j)

    def form_bond(self, atom_i, atom_j):
        # Mutual bond registration
        atom_i.bonds.append(atom_j)
        atom_j.bonds.append(atom_i)
        self.bonds.add(frozenset({atom_i.id, atom_j.id}))

        # Mass-weighted energy dissipation — simulates exothermic bonding
        v_rel = atom_i.vel - atom_j.vel
        total_mass = atom_i.mass + atom_j.mass
        delta_v = dissipation * v_rel  # fixed — was config.dissipation

        atom_i.vel -= (atom_j.mass / total_mass) * delta_v
        atom_j.vel += (atom_i.mass / total_mass) * delta_v

        # Record flash effect at bond midpoint
        midpoint = (atom_i.pos + atom_j.pos) / 2
        self.flash_events.append((midpoint.copy(), 0))

    def break_bond(self, atom_i, atom_j):
        # Remove each from the other's bond list
        atom_i.bonds.remove(atom_j)
        atom_j.bonds.remove(atom_i)

        # Remove pair from bond set — discard avoids KeyError if missing
        self.bonds.discard(frozenset({atom_i.id, atom_j.id}))


def detect_molecules(atoms):
    visited = set()
    molecule_counts = Counter()

    for atom in atoms:
        if atom.id in visited:
            continue

        # BFS to collect connected component (one molecule)
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

        # Count elements in this component
        counts = Counter(a.element for a in component)

        # Look up molecule name from config
        key = frozenset(counts.items())
        name = known_molecules.get(key, 'unknown')

        molecule_counts[name] += 1

    return molecule_counts
