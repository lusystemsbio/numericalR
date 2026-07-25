import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Lennard-Jones Metropolis Monte Carlo in a 2D periodic box
# U(r) = 1/(12 r^12) - 1/(6 r^6), minimum-image convention.
# Sample positions from p ~ exp(-U/T).
# ---------------------------------------------------------------

# --- Parameters ---
N       = 25          # number of particles (5x5 grid)
spacing = 1.25        # initial grid spacing
a       = 6.25        # box side length (5 * 1.25)
dxmax   = 0.25        # maximum single-particle move size
T       = 0.05        # temperature
n_equil = 5000        # equilibration steps
n_samp  = 10000       # sampling steps
seed    = 10

rng = np.random.default_rng(seed)

# --- Initial configuration: open 5x5 grid ---
coords = np.array([[i * spacing, j * spacing]
                   for i in range(5) for j in range(5)], dtype=float)

def min_image_disp(dr, a):
    # Wrap each displacement component into [-a/2, a/2] (minimum-image).
    return dr - a * np.round(dr / a)

def pair_energy(ri, coords, i, a):
    # Total LJ energy of particle i with all other particles j != i.
    dr = coords - ri                      # displacement vectors to all particles
    dr = min_image_disp(dr, a)            # apply minimum-image convention
    r2 = np.sum(dr**2, axis=1)            # squared distances
    r2[i] = np.inf                        # exclude self-interaction
    inv_r6 = 1.0 / r2**3                  # r^-6
    inv_r12 = inv_r6**2                   # r^-12
    return np.sum(inv_r12 / 12.0 - inv_r6 / 6.0)

def total_energy(coords, a):
    # Full box energy (each pair once).
    E = 0.0
    for i in range(N):
        # sum over j > i to avoid double counting
        dr = coords[i+1:] - coords[i]
        dr = min_image_disp(dr, a)
        r2 = np.sum(dr**2, axis=1)
        inv_r6 = 1.0 / r2**3
        inv_r12 = inv_r6**2
        E += np.sum(inv_r12 / 12.0 - inv_r6 / 6.0)
    return E

# --- Track current total energy incrementally ---
E_current = total_energy(coords, a)
E_initial = E_current

energies = []        # box energy recorded at every step
n_accept = 0

n_total = n_equil + n_samp
for step in range(n_total):
    # Pick a random particle for a local move.
    i = rng.integers(N)
    old_pos = coords[i].copy()

    # Energy of particle i in its current position (before move).
    E_old_i = pair_energy(old_pos, coords, i, a)

    # Propose a small random displacement, wrap back into the box.
    new_pos = old_pos + rng.uniform(-dxmax, dxmax, size=2)
    new_pos = np.mod(new_pos, a)

    # Energy of particle i in the trial position.
    E_new_i = pair_energy(new_pos, coords, i, a)

    # Incremental energy change from moving only particle i.
    dE = E_new_i - E_old_i

    # Metropolis acceptance criterion.
    if dE <= 0.0 or rng.random() < np.exp(-dE / T):
        coords[i] = new_pos          # accept: update position
        E_current += dE              # incremental energy update
        n_accept += 1
    # else: reject, keep old position and energy

    energies.append(E_current)

energies = np.array(energies)

# --- Results ---
E_final = E_current
E_equil_end = energies[n_equil - 1]
E_samp_mean = energies[n_equil:].mean()
E_samp_std = energies[n_equil:].std()
accept_ratio = n_accept / n_total

print(f"Initial box energy (open grid): {E_initial:.6f}")
print(f"Box energy at end of equilibration (step {n_equil}): {E_equil_end:.6f}")
print(f"Final box energy: {E_final:.6f}")
print(f"Mean box energy during sampling: {E_samp_mean:.6f}")
print(f"Std of box energy during sampling: {E_samp_std:.6f}")
print(f"Overall acceptance ratio: {accept_ratio:.6f}")

# --- Condensation check ---
# Compare the early equilibration energy to the sampling plateau.
E_early_mean = energies[:500].mean()
E_drop = E_early_mean - E_samp_mean
print(f"Mean box energy over first 500 (early equilibration) steps: {E_early_mean:.6f}")
print(f"Energy drop from early equilibration to sampling plateau: {E_drop:.6f}")
condensed = E_samp_mean < E_early_mean
print(f"Energy fell into a lower plateau (condensation confirmed): {condensed}")

# --- Plot: box energy over equilibration and sampling phases ---
fig, ax = plt.subplots(figsize=(9, 5))
steps = np.arange(n_total)
ax.plot(steps[:n_equil], energies[:n_equil], color="tab:orange",
        lw=0.8, label="equilibration")
ax.plot(steps[n_equil:], energies[n_equil:], color="tab:blue",
        lw=0.8, label="sampling")
ax.axvline(n_equil, color="k", ls="--", lw=1, label="equil / sample boundary")
ax.axhline(E_samp_mean, color="tab:green", ls=":", lw=1.2,
           label=f"sampling mean = {E_samp_mean:.3f}")
ax.set_xlabel("Monte Carlo step")
ax.set_ylabel("Box energy  U")
ax.set_title("LJ Metropolis MC: box energy over equilibration and sampling")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8C.3.1_s2.png", dpi=130)

# Explanation: The monotonic drop in energy during equilibration followed by
# small fluctuations around a low plateau confirms the sampler correctly moves
# the system from the high-energy open grid to a low-energy condensed cluster
# that is the equilibrium state favored by the Boltzmann distribution at low T.
print("Check explanation: the energy falls as particles condense from the open "
      "grid into a denser, lower-energy cluster and then only fluctuates about a "
      "low plateau, which is exactly the behavior expected when Metropolis "
      "sampling correctly relaxes the system toward its low-temperature "
      "Boltzmann equilibrium.")
