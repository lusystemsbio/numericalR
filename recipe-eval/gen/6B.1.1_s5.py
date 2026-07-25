import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
dx = 1.0        # step size in space
dt = 1.0        # step size in time
n_steps = 1000  # number of steps per walk
n_walks = 1000  # number of independent walkers
rng = np.random.default_rng(12)  # seed for reproducibility

# --- Method 1: explicit loop ---
# Build an (n_walks x (n_steps+1)) array; column 0 is the start at x=0.
walks_loop = np.zeros((n_walks, n_steps + 1))
for w in range(n_walks):              # loop over each independent walker
    x = 0.0                           # every walker starts at the origin
    for s in range(1, n_steps + 1):   # loop over time steps
        step = dx * (1 if rng.random() < 0.5 else -1)  # +dx or -dx with equal prob.
        x = x + step                  # x_next = x + dx*(+-1)
        walks_loop[w, s] = x          # record the new position

# --- Method 2: cumulative sum of random steps (vectorized) ---
# Draw all +/-1 steps at once, scale by dx, then accumulate along the time axis.
signs = rng.choice([-1.0, 1.0], size=(n_walks, n_steps))  # random +/-1 for every step
steps = dx * signs                                        # each move is +dx or -dx
positions = np.cumsum(steps, axis=1)                      # running total = position
walks_cumsum = np.hstack([np.zeros((n_walks, 1)), positions])  # prepend start at 0

# Time axis (0, dt, 2*dt, ...)
time = np.arange(n_steps + 1) * dt

# --- Plot several trajectories versus time ---
plt.figure(figsize=(9, 5))
for w in range(8):  # show 8 example walks to illustrate individual wandering
    plt.plot(time, walks_cumsum[w], lw=0.8, alpha=0.8)
plt.axhline(0, color="k", lw=0.6, ls="--")
plt.xlabel("time")
plt.ylabel("position x")
plt.title("Sample 1D random-walk trajectories")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6B.1.1_s5.png")

# --- Check: individual walks differ, but the ensemble spreads symmetrically ---
# Agreement between the two independent implementations (both use rng, so not identical
# paths, but both must obey the same statistics).
final_loop = walks_loop[:, -1]
final_cumsum = walks_cumsum[:, -1]

# Individual distinctness: how many of the sampled walks are pairwise identical.
unique_finals = len(np.unique(final_cumsum))

# Ensemble statistics at the final time.
mean_final = final_cumsum.mean()          # should be ~0 (symmetric about origin)
std_final = final_cumsum.std()            # spread; theory sqrt(n_steps)*dx
skew_final = (((final_cumsum - mean_final) / std_final) ** 3).mean()  # ~0 if symmetric
theory_std = np.sqrt(n_steps) * dx

print(f"Number of walks: {n_walks}")
print(f"Number of steps per walk: {n_steps}")
print(f"Distinct final positions among all walks: {unique_finals}")
print(f"Mean of final positions (loop method): {final_loop.mean():.6f}")
print(f"Mean of final positions (cumsum method): {mean_final:.6f}")
print(f"Std of final positions (cumsum method): {std_final:.6f}")
print(f"Theoretical std sqrt(N)*dx: {theory_std:.6f}")
print(f"Skewness of final positions: {skew_final:.6f}")
print(f"Min final position: {final_cumsum.min():.6f}")
print(f"Max final position: {final_cumsum.max():.6f}")

# One-sentence explanation:
# The mean staying near 0 with near-zero skewness while individual final positions
# take many distinct values confirms the result, because a mean/skew near zero means
# the cloud of walkers is balanced (symmetric) about the origin even though each walk
# follows its own unique, unpredictable path.
print("Explanation: a near-zero mean and skewness show the ensemble is symmetric about "
      "the origin, while the many distinct final positions show each walk wanders "
      "differently.")
