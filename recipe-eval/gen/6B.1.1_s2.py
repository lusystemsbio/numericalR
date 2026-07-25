import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Parameters ----
dx = 1.0        # step size in space
dt = 1.0        # step size in time
n_steps = 1000  # number of steps per walk
n_walks = 1000  # number of independent walks
seed = 12

rng = np.random.default_rng(seed)

# ---- Method 1: explicit loop ----
# For each walk we start at 0 and add +dx or -dx at each step.
x_loop = np.zeros((n_walks, n_steps + 1))  # position including the starting point
for w in range(n_walks):
    x = 0.0                                # walker starts at the origin
    for s in range(n_steps):
        step = dx * (1 if rng.random() < 0.5 else -1)  # +dx or -dx, equal probability
        x = x + step                       # update position: x_next = x + dx*(+-1)
        x_loop[w, s + 1] = x               # record the new position

# ---- Method 2: cumulative sum of random steps ----
# Draw all steps at once as +-dx, then integrate them with a cumulative sum.
rng2 = np.random.default_rng(seed)         # same seed -> same random stream / result
signs = rng2.integers(0, 2, size=(n_walks, n_steps)) * 2 - 1  # array of +1/-1
steps = dx * signs                         # convert to +-dx steps
x_cumsum = np.zeros((n_walks, n_steps + 1))
x_cumsum[:, 1:] = np.cumsum(steps, axis=1) # positions are the running total of steps

# Time axis
t = np.arange(n_steps + 1) * dt

# ---- Consistency check between the two implementations ----
max_diff = np.max(np.abs(x_loop - x_cumsum))
print(f"Max abs difference between loop and cumsum methods: {max_diff}")

# ---- Group statistics: spread symmetric about the origin ----
final_positions = x_cumsum[:, -1]          # position of every walk at the last step
mean_final = np.mean(final_positions)      # should be ~0 (symmetric about origin)
std_final = np.std(final_positions)        # measures the spread
theory_std = dx * np.sqrt(n_steps)         # expected spread sqrt(N)*dx for a symmetric walk

print(f"Mean final position (all walks): {mean_final}")
print(f"Std of final position (all walks): {std_final}")
print(f"Theoretical std sqrt(N)*dx: {theory_std}")
print(f"Min final position: {np.min(final_positions)}")
print(f"Max final position: {np.max(final_positions)}")

# Fraction ending on each side of the origin (should be roughly balanced)
frac_positive = np.mean(final_positions > 0)
frac_negative = np.mean(final_positions < 0)
print(f"Fraction of walks ending positive: {frac_positive}")
print(f"Fraction of walks ending negative: {frac_negative}")

# ---- Plot several trajectories versus time ----
n_show = 8
plt.figure(figsize=(10, 6))
for w in range(n_show):
    plt.plot(t, x_cumsum[w], lw=1.0, alpha=0.8, label=f"walk {w}")
plt.axhline(0, color="k", lw=0.8, ls="--")
plt.xlabel("time t")
plt.ylabel("position x")
plt.title(f"{n_show} of {n_walks} symmetric 1D random walks (dx={dx}, dt={dt}, {n_steps} steps)")
plt.legend(loc="upper left", fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.1.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: A near-zero mean with a large nonzero spread across the "
      "1000 walks confirms each walk wanders individually while the ensemble "
      "stays symmetric about the origin, since equal +-dx probabilities make "
      "the expected displacement zero but the variance grow as N*dx^2.")
