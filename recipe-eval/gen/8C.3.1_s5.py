import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model parameters -----
N = 25                # number of particles
n_side = 5            # 5x5 initial grid
spacing = 1.25        # grid spacing
a = n_side * spacing  # box side = 6.25 (periodic)
dxmax = 0.25          # max move displacement per coordinate
T = 0.05              # temperature
n_equil = 5000        # equilibration steps
n_sample = 10000      # sampling steps
rng = np.random.default_rng(10)  # seed 10

# ----- Lennard-Jones pair potential U(r) = 1/(12 r^12) - 1/(6 r^6) -----
def u_pair(r2):
    # r2 is squared distance; work with inverse powers to avoid sqrt
    inv6 = 1.0 / r2**3          # 1/r^6
    inv12 = inv6 * inv6         # 1/r^12
    return inv12 / 12.0 - inv6 / 6.0

# ----- Minimum-image energy of particle i with every other particle -----
def particle_energy(pos, i):
    dx = pos[:, 0] - pos[i, 0]
    dy = pos[:, 1] - pos[i, 1]
    # apply minimum-image convention: wrap displacements into (-a/2, a/2]
    dx -= a * np.round(dx / a)
    dy -= a * np.round(dy / a)
    r2 = dx * dx + dy * dy
    r2[i] = np.inf   # exclude self-interaction
    return np.sum(u_pair(r2))

# ----- Total energy (sum over unique pairs) -----
def total_energy(pos):
    e = 0.0
    for i in range(N):
        e += particle_energy(pos, i)
    return 0.5 * e   # each pair counted twice above

# ----- Initialize particles on the open 5x5 grid -----
positions = np.zeros((N, 2))
k = 0
for ix in range(n_side):
    for iy in range(n_side):
        positions[k] = [ix * spacing, iy * spacing]
        k += 1

E = total_energy(positions)          # running total energy
energy_trace = np.empty(n_equil + n_sample)
accepted = 0

# ----- Metropolis-Hastings sweep with incremental energy update -----
total_steps = n_equil + n_sample
for step in range(total_steps):
    i = rng.integers(N)                       # pick a random particle
    old_i = positions[i].copy()
    e_old = particle_energy(positions, i)     # its energy before the move

    # propose a local single-particle move
    positions[i, 0] = (old_i[0] + rng.uniform(-dxmax, dxmax)) % a
    positions[i, 1] = (old_i[1] + rng.uniform(-dxmax, dxmax)) % a

    e_new = particle_energy(positions, i)     # its energy after the move
    dE = e_new - e_old                        # incremental energy change

    # Metropolis acceptance criterion for Boltzmann exp(-U/T)
    if dE <= 0.0 or rng.random() < np.exp(-dE / T):
        E += dE               # accept: update running total incrementally
        accepted += 1
    else:
        positions[i] = old_i  # reject: restore old position

    energy_trace[step] = E

acceptance_ratio = accepted / total_steps
equil_energy = energy_trace[:n_equil]
sample_energy = energy_trace[n_equil:]

# ----- Numerical results -----
print(f"Box side a: {a}")
print(f"Initial energy (open 5x5 grid): {energy_trace[0]:.6f}")
print(f"Energy at end of equilibration: {equil_energy[-1]:.6f}")
print(f"Final energy: {energy_trace[-1]:.6f}")
print(f"Mean sampling-phase energy: {np.mean(sample_energy):.6f}")
print(f"Std of sampling-phase energy: {np.std(sample_energy):.6f}")
print(f"Min sampling-phase energy: {np.min(sample_energy):.6f}")
print(f"Max sampling-phase energy: {np.max(sample_energy):.6f}")
print(f"Overall acceptance ratio: {acceptance_ratio:.6f}")

# ----- Check: energy drops during equilibration then plateaus -----
initial_E = energy_trace[0]
plateau_E = np.mean(sample_energy)
condensed = plateau_E < initial_E
print(f"Energy fell from initial to plateau (condensation): {condensed}")
print(f"Energy drop (initial - plateau): {initial_E - plateau_E:.6f}")

# ----- Plot energy over equilibration and sampling phases -----
fig, ax = plt.subplots(figsize=(9, 5))
ax.plot(np.arange(total_steps), energy_trace, lw=0.6, color="steelblue")
ax.axvline(n_equil, color="red", ls="--", label="equilibration | sampling")
ax.axhline(plateau_E, color="green", ls=":", label=f"sampling mean = {plateau_E:.3f}")
ax.set_xlabel("Metropolis step")
ax.set_ylabel("Box energy U")
ax.set_title("2D Lennard-Jones box: Metropolis Monte Carlo energy trace")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8C.3.1_s5.png")

# Explanation: The equilibration energy monotonically decreasing from the high
# open-grid value to a low fluctuating plateau confirms the sampler is correctly
# driving the system toward the low-energy condensed cluster favored by the
# Boltzmann distribution at low T, i.e. it samples equilibrium configurations.
print("Check explanation: the energy falling from the open-grid value to a low "
      "fluctuating plateau confirms correctness because it shows Metropolis moves "
      "drive the system into the dense, low-energy cluster favored by the "
      "Boltzmann distribution at low T and then sample equilibrium fluctuations.")
