import pygame
import numpy as np
from config import *
from chemistry import *
from physics import *

def sim_to_screen(pos):
    return (int(pos[0] * scale), int(pos[1] * scale))

class Renderer:
    def __init__(self, screen):
        self.screen = screen
        self.large_font = pygame.font.SysFont('monospace', 18)
        self.small_font = pygame.font.SysFont('monospace', 13)

    def draw_background(self):
        # Draw simulation box
        pygame.draw.rect(self.screen, simulation_background_color, (0, 0, sim_width_px, height))
        # Draw palette panel
        pygame.draw.rect(self.screen, panel_background_color, (sim_width_px, 0, width - sim_width_px, height))

    def draw_bonds(self, atoms, bond_manager):
        # Build atom lookup dict for quick access by id
        atom_lookup = {atom.id: atom for atom in atoms}

        for pair in bond_manager.bonds:
            # Unpack pair ids and retrieve atoms
            id_i, id_j = tuple(pair)
            atom_i = atom_lookup[id_i]
            atom_j = atom_lookup[id_j]

            # Convert to screen coordinates as numpy arrays
            screen_i = np.array(sim_to_screen(atom_i.pos))
            screen_j = np.array(sim_to_screen(atom_j.pos))

            # Get bond order — default to 1 if not in bond_orders
            element_pair = frozenset({atom_i.element, atom_j.element})
            order = bond_orders.get(element_pair, 1)

            # Compute average color between the two atoms
            color_i = elements[atom_i.element]['color']
            color_j = elements[atom_j.element]['color']
            avg_color = tuple((np.array(color_i) + np.array(color_j)) // 2)

            # Compute perpendicular offset vector for multi-bond rendering
            bond_vec = screen_j - screen_i
            bond_len = np.linalg.norm(bond_vec)
            if bond_len > 0:
                perp = np.array([-bond_vec[1], bond_vec[0]]) / bond_len
            else:
                perp = np.array([0, 1])
            offset = perp * 3

            # Draw single bond
            if order == 1:
                pygame.draw.line(self.screen, avg_color,
                                tuple(screen_i.astype(int)),
                                tuple(screen_j.astype(int)), 2)
            # Draw double bond — two parallel lines
            elif order == 2:
                pygame.draw.line(self.screen, avg_color,
                                tuple((screen_i + offset).astype(int)),
                                tuple((screen_j + offset).astype(int)), 2)
                pygame.draw.line(self.screen, avg_color,
                                tuple((screen_i - offset).astype(int)),
                                tuple((screen_j - offset).astype(int)), 2)
            # Draw triple bond — center line plus two offset lines
            elif order == 3:
                pygame.draw.line(self.screen, avg_color,
                                tuple(screen_i.astype(int)),
                                tuple(screen_j.astype(int)), 2)
                pygame.draw.line(self.screen, avg_color,
                                tuple((screen_i + offset).astype(int)),
                                tuple((screen_j + offset).astype(int)), 2)
                pygame.draw.line(self.screen, avg_color,
                                tuple((screen_i - offset).astype(int)),
                                tuple((screen_j - offset).astype(int)), 2)

    def draw_atoms(self, atoms):
        for atom in atoms:
            # Convert position to screen coordinates
            screen_pos = sim_to_screen(atom.pos)
            color = elements[atom.element]['color']
            radius = elements[atom.element]['radius_pixel']

            # Draw ghost for held atoms
            if atom.held:
                ghost_surface = pygame.Surface((radius * 2 + 6, radius * 2 + 6), pygame.SRCALPHA)
                pygame.draw.circle(ghost_surface, (*color, 128), (radius + 3, radius + 3), radius)
                self.screen.blit(ghost_surface, (screen_pos[0] - radius - 3, screen_pos[1] - radius - 3))

            # Draw free atoms with outline and label
            else:
                # Filled circle — main atom body
                pygame.draw.circle(self.screen, color, screen_pos, radius)
                # White outline ring
                pygame.draw.circle(self.screen, (255, 255, 255), screen_pos, radius, 1)
                # Element label centered on atom
                text_surface = self.small_font.render(atom.element, True, (255, 255, 255))
                text_rect = text_surface.get_rect(center=screen_pos)
                self.screen.blit(text_surface, text_rect)

    def draw_flash_effects(self, bond_manager):
        for i in range(len(bond_manager.flash_events)):
            pos, age = bond_manager.flash_events[i]

            # Convert flash position to screen coordinates
            screen_pos = sim_to_screen(pos)

            # Expand radius and fade alpha over lifetime
            radius = int(age * 3)
            alpha = int(255 * (1 - age / 10))

            if radius > 0:
                size = radius * 2 + 4
                flash_surface = pygame.Surface((size, size), pygame.SRCALPHA)

                # Draw expanding ring
                pygame.draw.circle(
                    flash_surface,
                    (255, 255, 255, alpha),
                    (size // 2, size // 2),
                    radius,
                    2
                )

                # Blit centered on flash position
                offset = (screen_pos[0] - size // 2, screen_pos[1] - size // 2)
                self.screen.blit(flash_surface, offset)

            # Increment age for this flash event
            bond_manager.flash_events[i] = (pos, age + 1)

        # Remove expired flash events
        bond_manager.flash_events = [
            (pos, age) for (pos, age) in bond_manager.flash_events if age < 10
        ]

    def draw_panel(self, selected_element):
        # Compute horizontal center of panel
        panel_x_center = sim_width_px + (width - sim_width_px) // 2
        start_y = 80
        spacing = 120
        circle_radius = 25

        # Draw panel title
        title = self.large_font.render('Elements', True, (200, 200, 200))
        title_rect = title.get_rect(center=(panel_x_center, 30))
        self.screen.blit(title, title_rect)

        for i, element in enumerate(['H', 'O', 'N', 'He']):
            # Compute vertical position for this element button
            center_y = start_y + i * spacing
            color = elements[element]['color']

            # Draw element circle
            pygame.draw.circle(self.screen, color, (panel_x_center, center_y), circle_radius)

            # Draw selection ring if this element is selected
            if element == selected_element:
                pygame.draw.circle(self.screen, (255, 255, 255), (panel_x_center, center_y), circle_radius + 3, 2)

            # Draw element label below circle
            label = self.small_font.render(element, True, (255, 255, 255))
            label_rect = label.get_rect(center=(panel_x_center, center_y + circle_radius + 12))
            self.screen.blit(label, label_rect)

    def draw_hud(self, molecule_counts):
        # Draw HUD title
        title = self.large_font.render('Molecules', True, (200, 200, 200))
        self.screen.blit(title, (10, 10))

        y_offset = 35

        for name, count in molecule_counts.items():
            # Skip unrecognised molecules
            if name == 'unknown':
                continue

            # Draw molecule count line
            line = self.small_font.render(f"{name}: {count}", True, (200, 200, 200))
            self.screen.blit(line, (10, y_offset))
            y_offset += 18

    def draw(self, atoms, bond_manager, molecule_counts, selected_element):
        # Draw all layers in order
        self.draw_background()
        self.draw_bonds(atoms, bond_manager)
        self.draw_flash_effects(bond_manager)
        self.draw_atoms(atoms)
        self.draw_panel(selected_element)
        self.draw_hud(molecule_counts)

        # Flip display buffer
        pygame.display.flip()
