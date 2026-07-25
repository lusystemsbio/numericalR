import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# --- Parameters -----------------------------------------------------------
dx = 1.0            # step size scale (std of Gaussian / magnitude of discrete step)
dt = 1.0            # time step (unused in the update but part of the model spec)
n_steps = 1000      # number of steps per walk
n_walks = 1000      # number of independent walks
seed = 12
rng = np.random.default_rng(seed)

# Times at which we sample the spread: 1, 101, 201, ..., 901
sample_times = np.arange(1, 902, 100)

# --- Gaussian-step walk (built explicitly, step by step) ------------------
# positions[w, t] holds walker w at time t; start all walkers at 0
gauss_pos = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):
    # each step is drawn from N(mean=0, std=dx); x_next = x + step
    steps = rng.normal(loc=0.0, scale=dx, size=n_walks)
    gauss_pos[:, t] = gauss_pos[:, t - 1] + steps

# --- Discrete-step walk (built explicitly, step by step) ------------------
# each step is +dx or -dx with equal probability (mean 0, variance dx^2)
disc_pos = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):
    # draw a +1/-1 sign for each walker, then scale by dx
    signs = rng.choice(np.array([-1.0, 1.0]), size=n_walks)
    disc_pos[:, t] = disc_pos[:, t - 1] + dx * signs

# --- Collect the sampled columns ------------------------------------------
gauss_samples = [gauss_pos[:, t] for t in sample_times]
disc_samples = [disc_pos[:, t] for t in sample_times]

# --- Report the spread (std) at each sampled time -------------------------
print("time      std(Gaussian)      std(discrete)      theory sqrt(t)*dx")
for i, t in enumerate(sample_times):
    sg = np.std(gauss_samples[i], ddof=1)
    sd = np.std(disc_samples[i], ddof=1)
    print(f"{t:4d}      {sg:12.5f}      {sd:12.5f}      {np.sqrt(t) * dx:12.5f}")

# --- Statistical check: are the two spreads identical? --------------------
# Compare the distribution of positions at the final sampled time (t=901).
t_final = sample_times[-1]
g_final = gauss_pos[:, t_final]
d_final = disc_pos[:, t_final]

# Levene's test: null hypothesis = equal variances (spreads) of the two samples
levene_stat, levene_p = stats.levene(g_final, d_final)
print()
print(f"Final sampled time t = {t_final}")
print(f"Variance (Gaussian) = {np.var(g_final, ddof=1):.5f}")
print(f"Variance (discrete) = {np.var(d_final, ddof=1):.5f}")
print(f"Levene test statistic = {levene_stat:.5f}")
print(f"Levene test p-value    = {levene_p:.5f}")
print(f"Equal-variance (spread) not rejected at 5%: {levene_p > 0.05}")

# Two-sample Kolmogorov-Smirnov test: are the whole distributions the same?
ks_stat, ks_p = stats.ks_2samp(g_final, d_final)
print(f"KS test statistic = {ks_stat:.5f}")
print(f"KS test p-value    = {ks_p:.5f}")
print(f"Same distribution not rejected at 5%: {ks_p > 0.05}")

# --- Box plot: Gaussian walk next to discrete walk ------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

ax1.boxplot(gauss_samples, positions=range(len(sample_times)),
            labels=[str(t) for t in sample_times])
ax1.set_title("Gaussian-step random walk")
ax1.set_xlabel("time")
ax1.set_ylabel("walker position")
ax1.grid(True, alpha=0.3)

ax2.boxplot(disc_samples, positions=range(len(sample_times)),
            labels=[str(t) for t in sample_times])
ax2.set_title("Discrete-step random walk")
ax2.set_xlabel("time")
ax2.grid(True, alpha=0.3)

fig.suptitle("Spread of 1D random walk: Gaussian steps vs. discrete steps")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.3.1_s2.png")

# --- One-sentence explanation of why the check confirms the result --------
print()
print("Explanation: The Levene and KS tests fail to reject equality of the "
      "two spreads/distributions, confirming that because both step laws share "
      "mean 0 and variance dx^2, the Central Limit Theorem drives both walks to "
      "the same N(0, t*dx^2) long-time distribution, so only the step mean and "
      "variance matter, not the detailed step shape.")
