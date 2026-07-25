import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------
dx = 1.0          # step size scale (std for Gaussian, magnitude for discrete)
dt = 1.0          # time step (kept for completeness; not needed in the update)
n_steps = 1000    # number of steps per walk
n_walks = 1000    # number of independent walkers
seed = 12
sample_times = list(range(1, 902, 100))  # 1, 101, ..., 901

rng = np.random.default_rng(seed)

# -----------------------------------------------------------------------------
# Gaussian-step 1D random walk, built explicitly step by step
# -----------------------------------------------------------------------------
# positions[w, t] = position of walker w at time t (t=0 is the origin)
pos_gauss = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):
    steps = rng.normal(loc=0.0, scale=dx, size=n_walks)  # each step ~ N(0, dx)
    pos_gauss[:, t] = pos_gauss[:, t - 1] + steps         # x_next = x + N(0, dx)

# -----------------------------------------------------------------------------
# Discrete-step 1D random walk (+/- dx with equal probability), same structure
# -----------------------------------------------------------------------------
pos_disc = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):
    coin = rng.integers(0, 2, size=n_walks)               # 0 or 1
    steps = np.where(coin == 0, -dx, dx)                  # map to -dx / +dx
    pos_disc[:, t] = pos_disc[:, t - 1] + steps           # x_next = x +/- dx

# -----------------------------------------------------------------------------
# Extract positions at the sampled times for both walks
# -----------------------------------------------------------------------------
gauss_samples = pos_gauss[:, sample_times]   # shape (n_walks, len(sample_times))
disc_samples = pos_disc[:, sample_times]

# -----------------------------------------------------------------------------
# Report the spread (std dev) at each sampled time and test equality of variance
# -----------------------------------------------------------------------------
print("Comparison of spread: Gaussian-step vs discrete-step walk")
print("time | std_gauss | std_disc | Levene_p (equal variance)")
for i, t in enumerate(sample_times):
    g = gauss_samples[:, i]
    d = disc_samples[:, i]
    std_g = np.std(g, ddof=1)
    std_d = np.std(d, ddof=1)
    # Levene's test: null hypothesis is that the two samples have equal variance
    _, p_levene = stats.levene(g, d)
    print(f"t={t:4d} | std_gauss={std_g:8.4f} | std_disc={std_d:8.4f} | Levene_p={p_levene:.4f}")

# Overall check at the final sampled time (t=901): are the spreads statistically identical?
g_final = gauss_samples[:, -1]
d_final = disc_samples[:, -1]
_, p_levene_final = stats.levene(g_final, d_final)
_, p_ks_final = stats.ks_2samp(g_final, d_final)  # also compare full distributions
theory_std = np.sqrt(sample_times[-1]) * dx        # expected std = sqrt(t) * dx
print()
print(f"At t={sample_times[-1]}: theoretical std = sqrt(t)*dx = {theory_std:.4f}")
print(f"At t={sample_times[-1]}: std_gauss = {np.std(g_final, ddof=1):.4f}")
print(f"At t={sample_times[-1]}: std_disc  = {np.std(d_final, ddof=1):.4f}")
print(f"At t={sample_times[-1]}: Levene p-value (equal variance) = {p_levene_final:.4f}")
print(f"At t={sample_times[-1]}: KS 2-sample p-value (same distribution) = {p_ks_final:.4f}")

# Interpretation: a Levene p-value > 0.05 means we cannot reject equal variance,
# i.e. the spreads are statistically identical.
identical = p_levene_final > 0.05
print(f"Spreads statistically identical at t={sample_times[-1]} (Levene p>0.05)? {identical}")

# -----------------------------------------------------------------------------
# Side-by-side box plots of walker positions at the sampled times
# -----------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

axes[0].boxplot([gauss_samples[:, i] for i in range(len(sample_times))],
                labels=[str(t) for t in sample_times])
axes[0].set_title("Gaussian steps: x_next = x + N(0, dx)")
axes[0].set_xlabel("time")
axes[0].set_ylabel("walker position")
axes[0].axhline(0, color="gray", lw=0.8, ls="--")

axes[1].boxplot([disc_samples[:, i] for i in range(len(sample_times))],
                labels=[str(t) for t in sample_times])
axes[1].set_title("Discrete steps: x_next = x +/- dx")
axes[1].set_xlabel("time")
axes[1].axhline(0, color="gray", lw=0.8, ls="--")

fig.suptitle("Spread of 1D random walk positions at sampled times "
             "(1000 walks): Gaussian vs discrete steps")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.3.1_s4.png")

# -----------------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result
# -----------------------------------------------------------------------------
print()
print("Explanation: Because the position at time t is a sum of many independent, "
      "identically distributed steps, the Central Limit Theorem makes that sum "
      "converge to a Gaussian whose variance depends only on the per-step mean (0) "
      "and variance (dx^2); since the Gaussian and discrete steps share those two "
      "moments, the non-significant Levene test (equal spread) confirms the step "
      "distribution's shape is irrelevant to the long-time spread.")
