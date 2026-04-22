import pygame
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
        pygame.draw.rect(self.screen, simulation_background_color, (0, 0, sim_width_px, height))
        pygame.draw.rect(self.screen, panel_background_color, (sim_width_px, 0, width - sim_width_px, height))

    def draw_bonds(self, atoms, bond_manager):
        

        pass