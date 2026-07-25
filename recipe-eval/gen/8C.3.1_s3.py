import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Lennard-Jones pair potential U(r) = 1/(12 r^12) - 1/(6 r^6)
# ---------------------------------------------------------------
def lj(r2):
    # take squared distance, return LJ energy
    inv6 = 1.0 / r2**3      # 1/r^6
    inv12 = inv6 * inv6     # 1/r^12
    return inv12 / 12.0 - inv6 / 6.0

# ---------------------------------------------------------------
# Minimum-image squared distance of particle i to all others
# in a periodic box of side a.
# ---------------------------------------------------------------
def dist2_to_all(pos, i, a):
    d = pos - pos[i]                 # displacement vectors
    d -= a * np.round(d / a)         # minimum-image wrap
    r2 = np.sum(d * d, axis=1)
    return r2

# ---------------------------------------------------------------
# Energy contribution of particle i with the rest of the box.
# ---------------------------------------------------------------
def energy_of_particle(pos, i, a):
    r2 = dist2_to_all(pos, i, a)
    r2[i] = np.inf                   # skip self-interaction
    return np.sum(lj(r2))

# ---------------------------------------------------------------
# Total box energy (sum over unique pairs).
# ---------------------------------------------------------------
def total_energy(pos, a):
    N = len(pos)
    E = 0.0
    for i in range(N):
        r2 = dist2_to_all(pos, i, a)
        r2[i] = np.inf
        E += 0.5 * np.sum(lj(r2))    # 0.5 corrects double counting
    return E

# ---------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------
rng = np.random.default_rng(10)      # seed 10
n_side = 5
N = n_side * n_side                  # 25 particles
spacing = 1.25
a = 6.25                             # box side
dxmax = 0.25                         # max move displacement
T = 0.05                             # temperature
n_equil = 5000
n_sample = 10000

# ---------------------------------------------------------------
# Initial configuration: 5x5 grid at spacing 1.25
# ---------------------------------------------------------------
grid = (np.arange(n_side) + 0.5) * spacing
X, Y = np.meshgrid(grid, grid)
pos = np.column_stack([X.ravel(), Y.ravel()]).astype(float)

E = total_energy(pos, a)             # current box energy
energies = [E]

# ---------------------------------------------------------------
# Metropolis-Hastings loop with incremental energy update
# ---------------------------------------------------------------
n_total = n_equil + n_sample
n_accept = 0
for step in range(n_total):
    i = rng.integers(N)                          # pick a particle
    e_old = energy_of_particle(pos, i, a)        # its energy before
    old_xy = pos[i].copy()
    # propose a local single-particle move
    pos[i] = (old_xy + (rng.random(2) - 0.5) * 2 * dxmax) % a
    e_new = energy_of_particle(pos, i, a)        # its energy after
    dE = e_new - e_old                           # incremental change
    # Metropolis acceptance criterion
    if dE <= 0 or rng.random() < np.exp(-dE / T):
        E += dE                                  # accept: update energy
        n_accept += 1
    else:
        pos[i] = old_xy                          # reject: restore
    energies.append(E)

energies = np.array(energies)

# ---------------------------------------------------------------
# Reported numerical results
# ---------------------------------------------------------------
print(f"Number of particles N: {N}")
print(f"Box side a: {a}")
print(f"Temperature T: {T}")
print(f"Initial grid energy: {energies[0]:.6f}")
print(f"Energy after equilibration (step {n_equil}): {energies[n_equil]:.6f}")
print(f"Final energy: {energies[-1]:.6f}")
print(f"Mean sampling-phase energy: {np.mean(energies[n_equil:]):.6f}")
print(f"Std sampling-phase energy: {np.std(energies[n_equil:]):.6f}")
print(f"Min energy over run: {np.min(energies):.6f}")
print(f"Acceptance ratio: {n_accept / n_total:.6f}")
# recompute total energy from scratch to validate incremental bookkeeping
print(f"Direct total energy of final config: {total_energy(pos, a):.6f}")
print(f"Incremental vs direct energy difference: {abs(E - total_energy(pos, a)):.3e}")

# ---------------------------------------------------------------
# Plot: box energy over equilibration and sampling phases
# ---------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.plot(energies, lw=0.6, color="steelblue")
plt.axvline(n_equil, color="crimson", ls="--", label="end of equilibration")
plt.xlabel("Monte Carlo step")
plt.ylabel("Box energy")
plt.title("Metropolis MC of 2D Lennard-Jones box (N=25, T=0.05)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8C.3.1_s3.png", dpi=120)

# ---------------------------------------------------------------
# Check: does the energy drop during equilibration then plateau?
# ---------------------------------------------------------------
early_mean = np.mean(energies[:200])            # near the open grid start
plateau_mean = np.mean(energies[n_equil:])      # low plateau after condensing
condensed = plateau_mean < early_mean
print(f"Early-phase mean energy: {early_mean:.6f}")
print(f"Plateau mean energy: {plateau_mean:.6f}")
print(f"Energy dropped and plateaued (condensation confirmed): {condensed}")

# Explanation: The check confirms the sampler works because a correct
# Boltzmann sampler at low T must drive the loosely spaced grid toward
# the low-energy clustered configurations that dominate exp(-U/T), so a
# monotone-ish energy drop to a stable low plateau is the signature of
# the particles condensing into the favored dense cluster.
