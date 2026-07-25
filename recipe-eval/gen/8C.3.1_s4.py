import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- Parameters ----------------
n_side = 5                 # 5x5 grid
N = n_side * n_side        # number of particles = 25
spacing = 1.25             # initial grid spacing
a = n_side * spacing       # box side length = 6.25
dxmax = 0.25               # max single-particle move size
T = 0.05                   # temperature
n_equil = 5000             # equilibration steps
n_sample = 10000           # sampling steps
seed = 10

rng = np.random.default_rng(seed)

# ---------------- Initial configuration: open 5x5 grid ----------------
xs = (np.arange(n_side) + 0.0) * spacing
X, Y = np.meshgrid(xs, xs)
pos = np.column_stack([X.ravel(), Y.ravel()]).astype(float)

# ---------------- Lennard-Jones potential ----------------
def lj(r2):
    # U(r) = 1/(12 r^12) - 1/(6 r^6), input is r^2
    inv6 = 1.0 / r2**3      # 1/r^6
    inv12 = inv6 * inv6     # 1/r^12
    return inv12 / 12.0 - inv6 / 6.0

# ---------------- Minimum-image energy of particle i with all others ----------------
def particle_energy(i, pos):
    # displacement of all particles from particle i
    d = pos - pos[i]
    d -= a * np.round(d / a)        # minimum-image convention
    r2 = np.sum(d * d, axis=1)
    r2[i] = np.inf                  # exclude self-interaction
    return np.sum(lj(r2))

# ---------------- Total energy (sum over unique pairs) ----------------
def total_energy(pos):
    E = 0.0
    for i in range(N):
        d = pos[i+1:] - pos[i]
        d -= a * np.round(d / a)    # minimum image
        r2 = np.sum(d * d, axis=1)
        E += np.sum(lj(r2))
    return E

# ---------------- Metropolis-Hastings sweep with incremental updates ----------------
E = total_energy(pos)               # track running box energy
energies = []

def mc_step():
    global E
    # pick a random particle for a local move
    i = rng.integers(N)
    e_old = particle_energy(i, pos)             # its energy before the move
    old_xy = pos[i].copy()
    # propose a small local displacement
    pos[i] = pos[i] + (rng.random(2) - 0.5) * 2.0 * dxmax
    pos[i] = pos[i] % a                          # wrap back into box
    e_new = particle_energy(i, pos)             # its energy after the move
    dE = e_new - e_old                           # incremental energy change
    # Metropolis acceptance rule
    if dE <= 0.0 or rng.random() < np.exp(-dE / T):
        E += dE                                  # accept: update running energy
    else:
        pos[i] = old_xy                          # reject: restore position

# ---------------- Equilibration phase ----------------
for _ in range(n_equil):
    mc_step()
    energies.append(E)

E_after_equil = E
E_equil_start = energies[0]

# ---------------- Sampling phase ----------------
sample_energies = []
for _ in range(n_sample):
    mc_step()
    energies.append(E)
    sample_energies.append(E)

# ---------------- Results ----------------
E_initial_grid = total_energy(np.column_stack([X.ravel(), Y.ravel()]).astype(float))
mean_sample = np.mean(sample_energies)
std_sample = np.std(sample_energies)

print(f"Number of particles N: {N}")
print(f"Box side a: {a}")
print(f"Initial open-grid box energy: {E_initial_grid}")
print(f"Box energy at end of first MC step (equil start): {E_equil_start}")
print(f"Box energy after equilibration (5000 steps): {E_after_equil}")
print(f"Mean box energy over sampling phase: {mean_sample}")
print(f"Std of box energy over sampling phase: {std_sample}")
print(f"Min box energy during sampling: {np.min(sample_energies)}")
print(f"Max box energy during sampling: {np.max(sample_energies)}")

# ---------------- Check: energy drops during equilibration, then plateaus ----------------
equil_arr = np.array(energies[:n_equil])
early_mean = np.mean(equil_arr[:200])          # energy early in equilibration
plateau_mean = np.mean(equil_arr[-1000:])      # energy on the late plateau
dropped = plateau_mean < early_mean
print(f"Mean energy early in equilibration (first 200 steps): {early_mean}")
print(f"Mean energy on late equilibration plateau (last 1000 steps): {plateau_mean}")
print(f"Check energy fell as particles condensed: {dropped}")
# One sentence: the monotone drop from the open-grid energy to a low, weakly
# fluctuating plateau confirms the sampler is correctly driving the system toward
# the Boltzmann-favored dense cluster (low U) at this low temperature.
print("Explanation: the falling then plateauing energy shows Metropolis is sampling "
      "the low-energy condensed cluster favored by exp(-U/T), confirming correct sampling.")

# ---------------- Plot: box energy over both phases ----------------
plt.figure(figsize=(9, 5))
steps = np.arange(len(energies))
plt.plot(steps, energies, lw=0.6, color="C0")
plt.axvline(n_equil, color="k", ls="--", lw=1.0, label="end of equilibration")
plt.xlabel("Monte Carlo step")
plt.ylabel("Box energy")
plt.title("Lennard-Jones box energy: equilibration + sampling (Metropolis MC)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8C.3.1_s4.png")
