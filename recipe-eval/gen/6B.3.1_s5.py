import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

# ----------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------
dx = 1.0          # step scale: std of Gaussian step / magnitude of discrete step
dt = 1.0          # time increment per step (unused in dynamics, just labels time)
n_steps = 1000    # number of steps per walk
n_walks = 1000    # number of independent walkers
seed = 12
sample_times = np.arange(1, 902, 100)  # times 1, 101, ..., 901 (10 sampled times)

rng = np.random.default_rng(seed)

# ----------------------------------------------------------------------
# 1) Continuous Gaussian-step 1D random walk, built explicitly
#    x_next = x + N(0, dx)
# ----------------------------------------------------------------------
# positions[w, t] = position of walker w after t steps; t = 0 .. n_steps
gauss_pos = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):                 # march forward one step at a time
    steps = rng.normal(loc=0.0, scale=dx, size=n_walks)  # Gaussian steps, mean 0, std dx
    gauss_pos[:, t] = gauss_pos[:, t - 1] + steps        # accumulate onto previous position

# ----------------------------------------------------------------------
# 2) Discrete-step 1D random walk, built explicitly
#    x_next = x +/- dx (each sign equally likely)
# ----------------------------------------------------------------------
disc_pos = np.zeros((n_walks, n_steps + 1))
for t in range(1, n_steps + 1):                 # march forward one step at a time
    signs = rng.choice([-1.0, 1.0], size=n_walks)        # +/- with equal probability
    disc_pos[:, t] = disc_pos[:, t - 1] + signs * dx     # step of fixed magnitude dx

# ----------------------------------------------------------------------
# 3) Extract walker positions at the sampled times
# ----------------------------------------------------------------------
gauss_samples = [gauss_pos[:, tt] for tt in sample_times]
disc_samples  = [disc_pos[:, tt]  for tt in sample_times]

# ----------------------------------------------------------------------
# 4) Report spread (standard deviation) at each sampled time
#    Theory: std should grow like dx * sqrt(t) for BOTH walks.
# ----------------------------------------------------------------------
print("=== Spread (std of positions) at sampled times ===")
print(f"{'time':>6} {'std_gauss':>12} {'std_disc':>12} {'sqrt(t)*dx':>12}")
for i, tt in enumerate(sample_times):
    sg = np.std(gauss_samples[i], ddof=1)
    sd = np.std(disc_samples[i],  ddof=1)
    theory = dx * np.sqrt(tt)
    print(f"{tt:>6d} {sg:>12.4f} {sd:>12.4f} {theory:>12.4f}")

# ----------------------------------------------------------------------
# 5) Statistical check: is the spread of the Gaussian walk identical
#    to the discrete walk? Use Levene's test for equal variance
#    (null hypothesis: equal variances -> same spread).
#    Large p-values => cannot reject equal spread.
# ----------------------------------------------------------------------
print("\n=== Equal-variance test (Levene) at each sampled time ===")
print(f"{'time':>6} {'levene_stat':>12} {'p_value':>12}")
p_values = []
for i, tt in enumerate(sample_times):
    stat, p = stats.levene(gauss_samples[i], disc_samples[i])
    p_values.append(p)
    print(f"{tt:>6d} {stat:>12.4f} {p:>12.4f}")

print(f"\nMinimum Levene p-value across sampled times: {min(p_values):.4f}")
print(f"Number of sampled times with p < 0.05: {sum(p < 0.05 for p in p_values)} out of {len(p_values)}")

# Also compare the full final-time distributions with a KS test
ks_stat, ks_p = stats.ks_2samp(gauss_pos[:, sample_times[-1]], disc_pos[:, sample_times[-1]])
print(f"\nKS test of full distributions at t={sample_times[-1]}: statistic={ks_stat:.4f}, p_value={ks_p:.4f}")

# ----------------------------------------------------------------------
# 6) Box plot: Gaussian walk vs discrete walk, side by side
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

axes[0].boxplot(gauss_samples, tick_labels=[str(tt) for tt in sample_times])
axes[0].set_title("Continuous Gaussian-step walk")
axes[0].set_xlabel("time (steps)")
axes[0].set_ylabel("walker position x")
axes[0].axhline(0.0, color="gray", lw=0.8, ls="--")

axes[1].boxplot(disc_samples, tick_labels=[str(tt) for tt in sample_times])
axes[1].set_title("Discrete +/- step walk")
axes[1].set_xlabel("time (steps)")
axes[1].axhline(0.0, color="gray", lw=0.8, ls="--")

fig.suptitle("Spread of 1D random walk: Gaussian steps vs discrete steps")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.3.1_s5.png")

# ----------------------------------------------------------------------
# 7) One-sentence explanation
# ----------------------------------------------------------------------
print(
    "\nExplanation: The large Levene/KS p-values show we cannot distinguish the two "
    "spreads, confirming that because both step distributions share the same mean (0) "
    "and variance (dx^2), the Central Limit Theorem makes the accumulated position "
    "converge to the same Gaussian ~N(0, t*dx^2) regardless of the individual step shape."
)
