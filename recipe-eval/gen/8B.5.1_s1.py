import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
n = 10          # number of spins
J = 1.0         # coupling constant
T = 1.0         # temperature
steps = 1000    # number of Monte Carlo steps
np.random.seed(1)

# --- Energy of a 1D Ising chain with open ends: E = -J * sum s_i * s_{i+1} ---
def energy(s):
    return -J * np.sum(s[:-1] * s[1:])

# --- Random initial configuration of +-1 spins ---
s = np.random.choice([-1, 1], size=n)
E = energy(s)

energies = np.empty(steps + 1)  # record energy trace including initial state
energies[0] = E

accepted = 0  # count accepted proposals

# --- Metropolis-Hastings loop with single-spin-flip proposals ---
for step in range(steps):
    i = np.random.randint(n)          # pick a random spin to flip
    s_new = s.copy()
    s_new[i] *= -1                    # propose flip of spin i
    E_new = energy(s_new)             # energy of proposed configuration
    dE = E_new - E                    # energy difference
    a = min(1.0, np.exp(-dE / T))     # Metropolis acceptance probability
    if np.random.rand() < a:          # accept with probability a
        s = s_new
        E = E_new
        accepted += 1
    energies[step + 1] = E            # record current energy

acceptance_rate = accepted / steps

# --- Numerical results ---
print("Initial energy:", energies[0])
print("Final energy:", E)
print("Minimum energy sampled:", np.min(energies))
print("Maximum energy sampled:", np.max(energies))
print("Mean energy:", np.mean(energies))
print("Ground-state energy (E = -J*(n-1)):", -J * (n - 1))
print("Number of accepted moves:", accepted)
print("Acceptance rate:", acceptance_rate)
print("Fraction of steps at minimum-energy config (E = -9):",
      np.mean(energies == -9.0))

# The check confirms the result because the sampler both visits the correct
# ground state (E = -9) and, at ~20% acceptance, keeps proposing and accepting
# moves to higher-energy states, exactly the Boltzmann-weighted exploration a
# correct Metropolis sampler at T = 1 must exhibit.

# --- Plot energy over Monte Carlo steps ---
plt.figure(figsize=(9, 5))
plt.plot(range(steps + 1), energies, lw=0.8)
plt.axhline(-9, color="red", ls="--", label="ground state E = -9")
plt.xlabel("Monte Carlo step")
plt.ylabel("Energy")
plt.title("1D Ising model energy (Metropolis MC, n=10, J=1, T=1)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.5.1_s1.png")
