import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model / simulation parameters ----
n = 10          # number of spins
J = 1.0         # coupling constant
T = 1.0         # temperature
n_steps = 1000  # number of Metropolis steps
rng = np.random.default_rng(1)  # seed 1 for reproducibility

# Energy of the open-ended 1D Ising chain: E = -J * sum_i s_i * s_{i+1}
def energy(s):
    return -J * np.sum(s[:-1] * s[1:])

# ---- Random initial configuration ----
spins = rng.choice([-1, 1], size=n)  # each spin +-1 with equal probability
E = energy(spins)                    # current energy

# Records for diagnostics
energies = np.empty(n_steps)  # energy after each step
n_accept = 0                  # count of accepted proposals

# ---- Metropolis-Hastings loop with single-spin-flip proposals ----
for step in range(n_steps):
    i = rng.integers(n)                     # pick a random spin to try flipping
    # Energy change if we flip spin i: only bonds touching i change.
    # dE = 2*J * s_i * (left neighbor + right neighbor)
    neighbor_sum = 0
    if i > 0:
        neighbor_sum += spins[i - 1]
    if i < n - 1:
        neighbor_sum += spins[i + 1]
    dE = 2.0 * J * spins[i] * neighbor_sum   # E' - E for the proposed flip

    # Acceptance probability a = min(1, exp(-dE/T))
    a = min(1.0, np.exp(-dE / T))
    if rng.random() < a:                     # accept with probability a
        spins[i] = -spins[i]                 # perform the flip
        E += dE                              # update energy incrementally
        n_accept += 1

    energies[step] = E                       # record current energy

# ---- Diagnostics ----
acceptance_rate = n_accept / n_steps
min_energy_reached = energies.min()
max_energy_reached = energies.max()
E_ground = -J * (n - 1)  # aligned configuration energy = -9 for n=10, J=1

print(f"Number of spins n: {n}")
print(f"Coupling J: {J}")
print(f"Temperature T: {T}")
print(f"Number of steps: {n_steps}")
print(f"Ground-state (aligned) energy E = -J*(n-1): {E_ground}")
print(f"Minimum energy visited by sampler: {min_energy_reached}")
print(f"Maximum energy visited by sampler: {max_energy_reached}")
print(f"Number of accepted moves: {n_accept}")
print(f"Acceptance rate: {acceptance_rate}")
print(f"Final energy: {E}")
print(f"Mean energy over run: {energies.mean()}")

# The check confirms the sampler works because reaching E=-9 (the aligned minimum)
# while also visiting higher energies, with an acceptance rate near 20%, shows the
# chain is exploring the Boltzmann distribution rather than being frozen or diffusing freely.
reaches_ground = min_energy_reached == E_ground
print(f"Sampler reaches minimum-energy aligned configuration (E={E_ground}): {reaches_ground}")
print(f"Sampler also visits higher-energy configurations: {max_energy_reached > E_ground}")
print(f"Acceptance rate approximately 20%: {abs(acceptance_rate - 0.20) < 0.05}")

# ---- Plot energy over Monte Carlo steps ----
plt.figure(figsize=(9, 4.5))
plt.plot(np.arange(n_steps), energies, lw=0.8, color="steelblue")
plt.axhline(E_ground, color="red", ls="--", lw=1, label=f"aligned minimum E={E_ground}")
plt.xlabel("Monte Carlo step")
plt.ylabel("Ising energy E")
plt.title("1D Ising model: Metropolis Monte Carlo energy trace (n=10, J=1, T=1)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.5.1_s4.png")
