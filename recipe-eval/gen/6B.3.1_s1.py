import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ---- Parameters ----
dx = 1.0          # step scale: std of Gaussian step / magnitude of discrete step
dt = 1.0          # time increment per step (unused in position update, kept for clarity)
n_steps = 1000    # number of steps per walk
n_walks = 1000    # number of independent walkers
seed = 12
sample_times = list(range(1, 902, 100))  # times 1, 101, ..., 901

rng = np.random.default_rng(seed)

# ---- Gaussian-step walk (built explicitly, step by step) ----
# Each walker starts at 0; positions[w, t] = position of walker w after t steps.
gauss_pos = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):
    # draw one Gaussian step ~ N(0, dx) for every walker and add to previous position
    step = rng.normal(0.0, dx, size=n_walks)
    gauss_pos[:, t] = gauss_pos[:, t - 1] + step

# ---- Discrete-step walk (built explicitly, step by step) ----
# Each step is +dx or -dx with equal probability (mean 0, variance dx^2 -- same as Gaussian).
disc_pos = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):
    # draw +1/-1 for every walker, scale by dx, add to previous position
    step = rng.choice([-1.0, 1.0], size=n_walks) * dx
    disc_pos[:, t] = disc_pos[:, t - 1] + step

# ---- Spread (standard deviation) at each sampled time ----
gauss_std = np.array([gauss_pos[:, t].std(ddof=1) for t in sample_times])
disc_std  = np.array([disc_pos[:, t].std(ddof=1)  for t in sample_times])
theory_std = np.array([np.sqrt(t) * dx for t in sample_times])  # expected sqrt(t)*dx

for i, t in enumerate(sample_times):
    print(f"time {t:4d}: Gaussian std = {gauss_std[i]:.4f}, "
          f"discrete std = {disc_std[i]:.4f}, theory sqrt(t)*dx = {theory_std[i]:.4f}")

# ---- Statistical check that spreads are identical ----
# Compare the spread at the final sampled time with tests sensitive to variance/distribution.
t_last = sample_times[-1]
g_last = gauss_pos[:, t_last]
d_last = disc_pos[:, t_last]

# Levene's test: null hypothesis = equal variances (equal spread)
lev_stat, lev_p = stats.levene(g_last, d_last)
# F-test on the variances (ratio of variances under equal-variance null)
var_g = g_last.var(ddof=1)
var_d = d_last.var(ddof=1)
F = var_g / var_d
df1 = df2 = n_walks - 1
f_p = 2 * min(stats.f.cdf(F, df1, df2), 1 - stats.f.cdf(F, df1, df2))
# Two-sample Kolmogorov-Smirnov test on the full distributions
ks_stat, ks_p = stats.ks_2samp(g_last, d_last)

print(f"variance at time {t_last}: Gaussian = {var_g:.4f}, discrete = {var_d:.4f}")
print(f"Levene test (equal variance): statistic = {lev_stat:.4f}, p-value = {lev_p:.4f}")
print(f"F-test variance ratio F = {F:.4f}, p-value = {f_p:.4f}")
print(f"KS two-sample test: statistic = {ks_stat:.4f}, p-value = {ks_p:.4f}")
print(f"Equal spread (Levene p > 0.05)?  {lev_p > 0.05}")
print(f"Equal distribution (KS p > 0.05)? {ks_p > 0.05}")

# ---- Box plots side by side ----
fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
axes[0].boxplot([gauss_pos[:, t] for t in sample_times], positions=range(len(sample_times)),
                showfliers=False)
axes[0].set_xticks(range(len(sample_times)))
axes[0].set_xticklabels(sample_times)
axes[0].set_title("Gaussian steps: x_next = x + N(0, dx)")
axes[0].set_xlabel("time (step)")
axes[0].set_ylabel("walker position")

axes[1].boxplot([disc_pos[:, t] for t in sample_times], positions=range(len(sample_times)),
                showfliers=False)
axes[1].set_xticks(range(len(sample_times)))
axes[1].set_xticklabels(sample_times)
axes[1].set_title("Discrete steps: x_next = x +/- dx")
axes[1].set_xlabel("time (step)")

fig.suptitle("Spread of 1D random walk: Gaussian vs discrete steps")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.3.1_s1.png")

# ---- One-sentence explanation ----
print("Explanation: because both step distributions share the same mean (0) and variance (dx^2), "
      "the Central Limit Theorem makes the sum of many steps converge to the same Gaussian with "
      "variance t*dx^2, so the statistical tests failing to reject equal spread confirms that only "
      "the step mean and variance -- not the step distribution's shape -- govern long-time behavior.")
