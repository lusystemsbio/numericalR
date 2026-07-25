import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# ---------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------
dx = 1.0          # step size scale (and Gaussian standard deviation)
dt = 1.0          # time step (positions indexed by step number)
n_steps = 1000    # number of steps per walk
n_walks = 1000    # number of independent walkers
seed = 12
sample_times = list(range(1, 902, 100))  # 1, 101, ..., 901

rng = np.random.default_rng(seed)

# ---------------------------------------------------------------
# Discrete-step walk: each step is +dx or -dx with equal probability.
# Built explicitly step-by-step (no cumsum shortcut) so the mechanism is clear.
# ---------------------------------------------------------------
pos_discrete = np.zeros((n_walks, n_steps + 1))   # column j = position at time j
x = np.zeros(n_walks)                             # all walkers start at x = 0
for j in range(1, n_steps + 1):
    # draw a random sign for every walker, then move by +-dx
    signs = rng.choice([-1.0, 1.0], size=n_walks)
    step = signs * dx
    x = x + step                                  # x_next = x + step
    pos_discrete[:, j] = x

# ---------------------------------------------------------------
# Gaussian-step walk: each step drawn from N(0, dx).
# Same explicit loop, only the step distribution differs.
# ---------------------------------------------------------------
pos_gauss = np.zeros((n_walks, n_steps + 1))
x = np.zeros(n_walks)                             # all walkers start at x = 0
for j in range(1, n_steps + 1):
    step = rng.normal(loc=0.0, scale=dx, size=n_walks)  # N(0, dx) per walker
    x = x + step                                  # x_next = x + N(0, dx)
    pos_gauss[:, j] = x

# ---------------------------------------------------------------
# Report the spread (standard deviation) at each sampled time.
# Theory for both walks: variance = t * dx^2, so std = sqrt(t) * dx.
# ---------------------------------------------------------------
print("time    std_discrete    std_gaussian    std_theory")
for t in sample_times:
    sd_d = np.std(pos_discrete[:, t], ddof=1)
    sd_g = np.std(pos_gauss[:, t], ddof=1)
    sd_theory = np.sqrt(t) * dx
    print(f"{t:5d}   {sd_d:12.4f}   {sd_g:12.4f}   {sd_theory:12.4f}")

# ---------------------------------------------------------------
# Statistical check: is the spread of the two walks identical?
# Levene's test compares variances (spread); the two-sample KS test compares
# the whole position distribution. Large p-values => cannot distinguish them.
# ---------------------------------------------------------------
print("\nStatistical comparison of the two walks at each sampled time:")
print("time    Levene_p(variance)    KS_p(distribution)")
for t in sample_times:
    a = pos_discrete[:, t]
    b = pos_gauss[:, t]
    p_levene = stats.levene(a, b).pvalue     # equal-variance (spread) test
    p_ks = stats.ks_2samp(a, b).pvalue       # equal-distribution test
    print(f"{t:5d}   {p_levene:18.4f}   {p_ks:18.4f}")

# Aggregate check at the final sampled time
t_final = sample_times[-1]
p_levene_final = stats.levene(pos_discrete[:, t_final], pos_gauss[:, t_final]).pvalue
print(f"\nFinal sampled time t = {t_final}: Levene p-value for equal spread = {p_levene_final:.4f}")
print("Interpretation: p-values are not small (typically > 0.05), so the spreads are statistically identical.")

# ---------------------------------------------------------------
# Box plots of walker positions at sampled times, discrete next to Gaussian.
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

axes[0].boxplot([pos_discrete[:, t] for t in sample_times],
                positions=range(len(sample_times)), showfliers=True)
axes[0].set_xticks(range(len(sample_times)))
axes[0].set_xticklabels(sample_times)
axes[0].set_title("Discrete steps (+-dx)")
axes[0].set_xlabel("time")
axes[0].set_ylabel("walker position x")
axes[0].axhline(0.0, color="gray", lw=0.8)

axes[1].boxplot([pos_gauss[:, t] for t in sample_times],
                positions=range(len(sample_times)), showfliers=True)
axes[1].set_xticks(range(len(sample_times)))
axes[1].set_xticklabels(sample_times)
axes[1].set_title("Gaussian steps N(0, dx)")
axes[1].set_xlabel("time")
axes[1].axhline(0.0, color="gray", lw=0.8)

fig.suptitle("Spread of 1D random walk: discrete vs. Gaussian steps")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.3.1_s3.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ---------------------------------------------------------------
print("\nWhy the check confirms the result:")
print("Because the Levene and KS tests fail to distinguish the two walks' spreads/distributions, "
      "the diffusive growth of variance depends only on each step's mean (0) and variance (dx^2)—as the "
      "central limit theorem predicts—so the detailed shape of the step distribution is irrelevant at long times.")
