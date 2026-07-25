import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Parameters --------------------------------------------------------
dx = 1.0        # spatial step size
dt = 1.0        # time step (used to build the time axis)
n_steps = 1000  # number of steps per walk
n_walks = 1000  # number of independent walkers
np.random.seed(12)

# ---- Method 1: explicit loop ------------------------------------------
# Each walk starts at 0. At every step we draw a random +1 or -1
# and add dx*(+-1) to the current position, recording the full path.
paths_loop = np.zeros((n_walks, n_steps + 1))  # column 0 is the start (x=0)
for w in range(n_walks):
    x = 0.0                                    # walker starts at the origin
    for s in range(1, n_steps + 1):
        step = dx * (1 if np.random.random() < 0.5 else -1)  # +dx or -dx
        x = x + step                           # x_next = x + dx*(+-1)
        paths_loop[w, s] = x                   # store the new position

# ---- Method 2: cumulative sum of random steps -------------------------
# Same model, vectorized: draw all +-1 steps at once, scale by dx,
# then take the running (cumulative) sum along time to get positions.
np.random.seed(12)                             # reseed so both methods match
steps = dx * np.where(np.random.random((n_walks, n_steps)) < 0.5, 1.0, -1.0)
paths_cumsum = np.zeros((n_walks, n_steps + 1))
paths_cumsum[:, 1:] = np.cumsum(steps, axis=1)  # accumulate steps over time

# ---- Cross-check that both implementations agree ----------------------
max_diff = np.max(np.abs(paths_loop - paths_cumsum))
print(f"Max abs difference between loop and cumsum methods: {max_diff}")

# ---- Time axis --------------------------------------------------------
time = np.arange(n_steps + 1) * dt

# ---- Statistics: individual variation vs. group symmetry --------------
final_positions = paths_cumsum[:, -1]          # position of each walk at the end
mean_final = np.mean(final_positions)          # should be ~0 (symmetric spread)
std_final = np.std(final_positions)            # spread of the ensemble
theory_std = dx * np.sqrt(n_steps)             # expected sqrt(N)*dx growth
median_final = np.median(final_positions)
frac_positive = np.mean(final_positions > 0)   # should be ~0.5 by symmetry

# Two specific walks to show they wander differently
diff_two_walks = np.max(np.abs(paths_cumsum[0] - paths_cumsum[1]))

print(f"Mean of final positions: {mean_final}")
print(f"Median of final positions: {median_final}")
print(f"Std of final positions (simulated): {std_final}")
print(f"Std of final positions (theory dx*sqrt(N)): {theory_std}")
print(f"Fraction of walks ending positive: {frac_positive}")
print(f"Max separation between walk 0 and walk 1: {diff_two_walks}")
print(f"Min final position: {np.min(final_positions)}")
print(f"Max final position: {np.max(final_positions)}")

# ---- Plot several trajectories ----------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
for w in range(8):                             # plot a handful of walks
    ax.plot(time, paths_cumsum[w], lw=1.0, alpha=0.8)
ax.axhline(0, color="k", lw=0.8, ls="--")       # reference: the origin
ax.set_xlabel("time (steps * dt)")
ax.set_ylabel("position x")
ax.set_title("Several 1D random-walk trajectories vs. time")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.1.1_s4.png")

# Explanation: The check confirms the result because the near-zero mean/median
# and ~50% positive fraction show the ensemble stays symmetric about the origin,
# while the large spread (std ~= dx*sqrt(N)) and the big separation between two
# individual walks show each walker follows its own distinct, unpredictable path.
print("Check: mean~0 and frac_positive~0.5 => symmetric spread about 0, "
      "while nonzero std and walk-to-walk separation => each walk wanders differently.")
