import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
dx = 1.0          # step size in space
dt = 1.0          # time increment per step
n_steps = 1000    # number of steps per walk
n_walks = 1000    # number of independent walkers
rng = np.random.default_rng(12)  # seeded RNG for reproducibility

# --- Method 1: explicit loop ---
# We build an (n_walks, n_steps+1) array, position 0 at t=0 for every walk.
pos_loop = np.zeros((n_walks, n_steps + 1))
# Draw all random +-1 steps up front: choice of {-1, +1} with equal probability.
steps = rng.choice([-1.0, 1.0], size=(n_walks, n_steps))
for w in range(n_walks):          # loop over each independent walker
    x = 0.0                       # every walker starts at the origin
    for s in range(n_steps):      # loop over time steps
        x = x + dx * steps[w, s]  # x_next = x + dx*(+-1)
        pos_loop[w, s + 1] = x    # record position after this step

# --- Method 2: cumulative sum of the same random steps ---
# cumsum along the time axis reproduces the loop result in one vectorized op.
pos_cumsum = np.zeros((n_walks, n_steps + 1))
pos_cumsum[:, 1:] = np.cumsum(dx * steps, axis=1)

# Confirm the two implementations agree exactly.
max_diff = np.max(np.abs(pos_loop - pos_cumsum))
print(f"Max abs difference between loop and cumsum implementations: {max_diff}")

# Use the cumsum result as the ensemble from here on.
positions = pos_cumsum
time = np.arange(n_steps + 1) * dt  # time axis

# --- Ensemble statistics at the final time ---
final = positions[:, -1]
mean_final = np.mean(final)                 # should be near 0 (symmetric spread)
std_final = np.std(final)                   # spread of the ensemble
theory_std = dx * np.sqrt(n_steps)          # expected sqrt(N) diffusive spread
skew_final = np.mean((final - mean_final)**3) / std_final**3  # ~0 if symmetric

# How different are individual walks? Compare spread across walks vs a single walk.
pairwise_example = np.abs(positions[0, -1] - positions[1, -1])

print(f"Number of walks: {n_walks}")
print(f"Number of steps: {n_steps}")
print(f"Final-position ensemble mean: {mean_final}")
print(f"Final-position ensemble std: {std_final}")
print(f"Theoretical std dx*sqrt(N): {theory_std}")
print(f"Final-position ensemble skewness: {skew_final}")
print(f"Final position of walk 0: {positions[0, -1]}")
print(f"Final position of walk 1: {positions[1, -1]}")
print(f"|walk0_final - walk1_final|: {pairwise_example}")
print(f"Fraction of walks ending x > 0: {np.mean(final > 0)}")
print(f"Fraction of walks ending x < 0: {np.mean(final < 0)}")

# --- Plot several trajectories versus time ---
plt.figure(figsize=(10, 6))
for w in range(8):  # plot a handful of individual trajectories
    plt.plot(time, positions[w], lw=1.0, alpha=0.8, label=f"walk {w}")
plt.axhline(0, color="k", lw=0.8, ls="--")
plt.xlabel("time")
plt.ylabel("position x")
plt.title("Sample 1D random-walk trajectories (dx=1, dt=1)")
plt.legend(loc="upper left", fontsize=8, ncol=2)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.1.1_s1.png")

# Explanation: The near-zero ensemble mean and skewness show the group is centered
# and symmetric about the origin, while each walk's distinct final position (and the
# large pairwise difference) shows the walkers wander independently -- confirming that
# individually random paths collectively spread symmetrically as expected for a
# discrete symmetric 1D random walk.
print("Check: a near-zero ensemble mean/skew with distinct individual final positions confirms that the walks wander differently yet spread symmetrically about the origin.")
