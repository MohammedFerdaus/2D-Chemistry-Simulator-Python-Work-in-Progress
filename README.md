# 2D Chemistry Simulator (Work In Progess)

A real-time 2D molecular dynamics simulator built in Python using pygame and numpy. Atoms of H, O, N, and He are drag-and-placed into a simulation box where Lennard-Jones forces, harmonic bond springs, Langevin thermostat noise, and angle constraints govern their motion. Bonds form and break dynamically, molecules are detected in real time via graph traversal, and the HUD displays all recognised chemical species as they appear.

---

## Repository Structure

```
chemistry_sim/
├── config.py       — all constants: display, physics, elements, bonding, molecules
├── physics.py      — Atom class, LJ force, spring force, Langevin thermostat, angle force, Velocity Verlet integrator
├── chemistry.py    — BondManager, bond formation/breaking, molecule detection via BFS
├── renderer.py     — all pygame drawing: atoms, bonds, flash effects, palette panel, HUD
└── main.py         — game loop, drag-and-place state machine, event handling
```

---

## Demo

### Drag and place
Place atoms into the sim box by selecting an element from the palette, clicking and dragging — release velocity sets the initial speed and direction.

![drag_and_place](gifs/drag_and_place.gif)

---

### Bond formation
When two atoms with open valence slots come within bonding range, a harmonic spring forms between them. A flash ring marks the moment of bonding and the molecule HUD updates in real time.

![bonding](gifs/bonding.gif)

---

### Repulsion
Unbound atoms interact via the Lennard-Jones potential — they repel strongly at short range and attract weakly at medium range, producing realistic collision behaviour.

![repulsion](gifs/repulsion.gif)

---

## Physics

### Lennard-Jones potential

Non-bonded atom pairs interact via the LJ potential, which combines a short-range repulsive wall with a medium-range attractive well:

$$V(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]$$

The scalar force derived from this potential is:

$$F(r) = \frac{24\varepsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]$$

Mixed-pair sigma is computed via the Lorentz-Berthelot combining rule:

$$\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2}$$

Forces are capped at `force_cap` and only computed within a cutoff radius $r_{cut} = 2.5\,\sigma_{ij}$.

---

### Harmonic bond spring

Bonded atom pairs are held together by a harmonic spring with equilibrium length $r_0$ set per element pair:

$$\vec{F}_{spring} = k(r - r_0)\hat{u}$$

If the bond is stretched beyond the break threshold $r_{break} = 2.5\,r_0$, the bond is removed and the pair returns to LJ interaction.

---

### Langevin thermostat

Each free atom receives a drag force and a stochastic kick every step, implementing the fluctuation-dissipation theorem to maintain temperature $T$:

$$\vec{F}_{drag} = -\gamma m \vec{v}$$

$$\vec{F}_{rand} = \sqrt{\frac{2\gamma m k_B T}{\Delta t}}\,\vec{\eta}, \quad \vec{\eta} \sim \mathcal{N}(0,1)^2$$

The balance between drag and noise keeps the system at a controlled thermal energy. Raising $T$ increases atomic speed; lowering it causes atoms to settle.

---

### Angle constraints

To enforce molecular geometry, a restoring torque is applied to every central atom with two or more bonds. For a central atom $B$ with bonded neighbours $A_1$ and $A_2$:

$$\cos\theta = \frac{\vec{r}_1 \cdot \vec{r}_2}{|\vec{r}_1||\vec{r}_2|}$$

$$\tau = -k_\theta(\theta - \theta_0)$$

Forces are applied perpendicular to each bond vector, pushing the angle toward the target $\theta_0$:

| Molecule | Central atom | $\theta_0$ |
|----------|-------------|-----------|
| H₂O | O | 104.5° |
| NH₃ | N | 107.0° |

---

### Velocity Verlet integrator

The integrator is second-order accurate and conserves energy far better than Euler. Each timestep runs in the following sequence:

$$\vec{r}(t+\Delta t) = \vec{r}(t) + \vec{v}(t)\Delta t + \frac{1}{2}\vec{a}(t)\Delta t^2$$

$$\vec{v}(t+\Delta t) = \vec{v}(t) + \frac{1}{2}\left[\vec{a}(t) + \vec{a}(t+\Delta t)\right]\Delta t$$

Forces are recomputed between the two half-steps so the velocity update uses the average of old and new accelerations.

---

## Elements

| Element | Max bonds | Role |
|---------|-----------|------|
| H | 1 | Forms bonds with O and N |
| O | 2 | Central atom in H₂O, H₂O₂ |
| N | 3 | Central atom in NH₃, N₂ |
| He | 0 | Noble gas — never bonds |

He is included specifically to demonstrate inert behaviour — it collides elastically with everything but never forms a bond regardless of proximity.

---

## Recognised Molecules

| Formula | Name |
|---------|------|
| H₂ | Hydrogen gas |
| O₂ | Oxygen gas |
| N₂ | Nitrogen gas |
| H₂O | Water |
| H₂O₂ | Hydrogen peroxide |
| OH· | Hydroxyl radical |
| HO₂· | Hydroperoxyl radical |
| NH₃ | Ammonia |
| N₂H₄ | Hydrazine |
| NH₂· | Amino radical |
| NO· | Nitric oxide |
| NO₂· | Nitrogen dioxide |
| N₂O | Nitrous oxide |
| N₂O₃ | Dinitrogen trioxide |
| HNO₃ | Nitric acid |
| HNO₂ | Nitrous acid |
| NH₂OH | Hydroxylamine |
| HNO | Nitroxyl |

---

## Stack

| Area | Library |
|------|---------|
| Rendering | pygame 2.6.1 |
| Numerics | numpy |
| Core | Python 3.13 standard library |

---

## How to Run

**Requirements:** Python 3.10+, pygame, numpy

Install dependencies:

```
pip install pygame numpy
```

Run the simulator:

```
python main.py
```

A 1100×700 window opens with the sim box on the left and the element palette on the right. Select an element from the palette, click and drag into the sim box to place atoms. Space to pause, Escape to quit.

---

## Controls

| Input | Action |
|-------|--------|
| Click element in palette | Select element |
| Click and drag in sim box | Place atom with velocity |
| Space | Pause / unpause |
| Escape | Quit |

---

## Tuning

All physics constants live in `config.py` and are the intended way to experiment:

| Constant | Effect |
|----------|--------|
| `temperature` | Controls thermal speed of atoms |
| `gamma` | Drag coefficient — higher values damp motion faster |
| `spring_k` | Bond stiffness — lower values reduce oscillation on bonding |
| `angle_spring_k` | Strength of geometry enforcement |
| `dt` | Timestep — reduce if simulation explodes |
| `epsilon` | LJ well depth — higher values strengthen repulsion |

---

## Notes

Built and tested on Python 3.13.2, pygame 2.6.1, Windows 10. Written entirely from scratch with no starter code. All physics implemented manually — no physics engine used.
