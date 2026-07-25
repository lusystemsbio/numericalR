import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Parameters ----
dx = 1.0          # step size in space
dt = 1.0          # time step
n_steps = 1000    # number of steps per walk
n_walks = 1000    # number of independent walks
seed = 12

rng = np.random.default_rng(seed)

# ---- Method 1: explicit loop ----
# Each walk starts at 0; at every step add +dx or -dx with equal probability.
positions_loop = np.zeros((n_walks, n_steps + 1))  # column 0 is the start (x=0)
for w in range(n_walks):            # loop over independent walkers
    x = 0.0                         # every walker starts at the origin
    for s in range(1, n_steps + 1): # loop over time steps
        step = dx if rng.random() < 0.5 else -dx  # +dx or -dx, 50/50
        x = x + step                # x_next = x + dx*(+/-1)
        positions_loop[w, s] = x    # record the new position

# ---- Method 2: cumulative sum of random steps (same idea, vectorized) ----
# Draw +/-1 for every (walk, step), scale by dx, then integrate with cumsum.
signs = rng.choice([-1.0, 1.0], size=(n_walks, n_steps))  # random +/-1
steps = dx * signs                                        # actual displacements
positions_cumsum = np.zeros((n_walks, n_steps + 1))       # column 0 stays at 0
positions_cumsum[:, 1:] = np.cumsum(steps, axis=1)        # running total = position

# Time axis
time = np.arange(n_steps + 1) * dt

# ---- Check A: each walk wanders differently ----
# Compare a few individual final positions; they should not all be the same.
final_positions = positions_cumsum[:, -1]
print("Final position of walk 0:", final_positions[0])
print("Final position of walk 1:", final_positions[1])
print("Final position of walk 2:", final_positions[2])
print("Number of distinct final positions among all walks:", len(np.unique(final_positions)))

# ---- Check B: as a group they spread symmetrically about the origin ----
# Ensemble mean should be ~0 (symmetric spreading); variance grows like n_steps*dx^2.
ensemble_mean_final = np.mean(final_positions)
ensemble_std_final = np.std(final_positions)
theoretical_std_final = np.sqrt(n_steps) * dx  # sqrt(N)*dx for a symmetric walk
frac_positive = np.mean(final_positions > 0)   # should be ~0.5 by symmetry
frac_negative = np.mean(final_positions < 0)

print("Ensemble mean of final position:", ensemble_mean_final)
print("Ensemble std of final position (simulated):", ensemble_std_final)
print("Ensemble std of final position (theory sqrt(N)*dx):", theoretical_std_final)
print("Fraction of walks ending positive:", frac_positive)
print("Fraction of walks ending negative:", frac_negative)

# ---- Consistency check: loop vs cumsum give statistically the same spread ----
print("Mean final position (loop method):", np.mean(positions_loop[:, -1]))
print("Std final position (loop method):", np.std(positions_loop[:, -1]))

# ---- Plot: several trajectories versus time ----
plt.figure(figsize=(9, 5))
for w in range(8):  # show 8 example trajectories
    plt.plot(time, positions_cumsum[w], lw=1.0, alpha=0.8)
plt.axhline(0.0, color="k", lw=0.8, ls="--")
plt.xlabel("time (t)")
plt.ylabel("position (x)")
plt.title("Several 1D symmetric random-walk trajectories")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.1.1_s3.png")

# Explanation:
# The near-zero ensemble mean with a roughly 50/50 split of positive/negative
# endpoints confirms symmetric spreading, because the +dx and -dx steps are
# equally likely, so individual walks diverge while their average cancels to zero.
print("Explanation: because +dx and -dx are equally likely, individual walks diverge yet the ensemble mean stays ~0 with ~50/50 positive/negative endpoints, confirming symmetric spreading about the origin.")
