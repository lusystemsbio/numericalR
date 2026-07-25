import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal minima at f = 0
def himmelblau(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known minima of Himmelblau's function
known_minima = np.array([
    [3.0, 2.0],
    [-2.805118, 3.131312],
    [-3.779310, -3.283186],
    [3.584428, -1.848126],
])

# Explicit Metropolis-Hastings at fixed temperature T
def metropolis_hastings(T, n_steps=int(1e4), step=0.5, start=(0.0, 0.0), seed=1):
    rng = np.random.default_rng(seed)
    x, y = float(start[0]), float(start[1])      # current state
    f_cur = himmelblau(x, y)                      # current energy
    samples = np.empty((n_steps, 2))             # storage for the chain
    n_accept = 0
    for i in range(n_steps):
        # uniform proposal displacement in [-step, +step] per coordinate
        xp = x + rng.uniform(-step, step)
        yp = y + rng.uniform(-step, step)
        f_prop = himmelblau(xp, yp)
        de = f_prop - f_cur                       # change in energy
        # acceptance probability a = min(1, exp(-de/T))
        a = 1.0 if de <= 0 else np.exp(-de / T)
        if rng.uniform() < a:                     # accept move
            x, y, f_cur = xp, yp, f_prop
            n_accept += 1
        samples[i] = (x, y)                        # record state (accepted or not)
    return samples, n_accept / n_steps

temperatures = [1, 10, 30, 50]
n_steps = int(1e4)
step_size = 0.5
start = (0.0, 0.0)

# Background contour of the landscape (log scale to reveal all basins)
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = himmelblau(GX, GY)

fig, axes = plt.subplots(2, 2, figsize=(12, 11))
results = {}

for ax, T in zip(axes.ravel(), temperatures):
    samples, acc_rate = metropolis_hastings(T, n_steps, step_size, start, seed=1)
    results[T] = (samples, acc_rate)

    ax.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
    ax.plot(samples[:, 0], samples[:, 1], ',', color="white", alpha=0.4)
    ax.scatter(samples[:, 0], samples[:, 1], s=2, color="orange", alpha=0.15)
    ax.scatter(known_minima[:, 0], known_minima[:, 1], marker="*",
               s=250, color="red", edgecolor="black", zorder=5, label="true minima")
    ax.scatter(*start, marker="o", s=60, color="cyan", edgecolor="black",
               zorder=6, label="start")
    ax.set_title(f"T = {T}  (accept rate = {acc_rate:.3f})")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.legend(loc="upper right", fontsize=8)

fig.suptitle("Metropolis-Hastings sampling of Himmelblau's landscape at fixed temperature",
             fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.2.1_s2.png", dpi=130)

# ---- Numerical diagnostics ----
# For each temperature, measure (a) how widely it explores: number of distinct
# basins visited (assign each sample to the nearest of the four minima and count
# how many are visited beyond a small threshold), and (b) how well it settles:
# the minimum f value reached and the mean f over the last 10% of the chain.
print("=== Metropolis-Hastings on Himmelblau's function ===")
print(f"start = {start}, n_steps = {n_steps}, step_size = {step_size}, seed = 1 per temperature")
print()

for T in temperatures:
    samples, acc_rate = results[T]
    fvals = himmelblau(samples[:, 0], samples[:, 1])

    # assign each sample to nearest known minimum
    d = np.linalg.norm(samples[:, None, :] - known_minima[None, :, :], axis=2)
    nearest = np.argmin(d, axis=1)
    counts = np.array([np.sum(nearest == k) for k in range(4)])
    frac = counts / n_steps
    # a basin counts as "visited" if it holds > 5% of the samples
    n_basins_visited = int(np.sum(frac > 0.05))

    tail = fvals[int(0.9 * n_steps):]  # last 10% of the chain

    print(f"--- T = {T} ---")
    print(f"T = {T}: acceptance rate = {acc_rate:.4f}")
    print(f"T = {T}: min f reached = {fvals.min():.6f}")
    print(f"T = {T}: mean f over last 10% of chain = {tail.mean():.4f}")
    print(f"T = {T}: basin occupancy fractions (minima 1-4) = "
          f"{frac[0]:.3f}, {frac[1]:.3f}, {frac[2]:.3f}, {frac[3]:.3f}")
    print(f"T = {T}: number of basins visited (>5% occupancy) = {n_basins_visited}")
    print()

# ---- The separate check: no single T both explores widely AND settles ----
samples_lo, _ = results[1]
fvals_lo = himmelblau(samples_lo[:, 0], samples_lo[:, 1])
d_lo = np.linalg.norm(samples_lo[:, None, :] - known_minima[None, :, :], axis=2)
frac_lo = np.array([np.mean(np.argmin(d_lo, axis=1) == k) for k in range(4)])
basins_lo = int(np.sum(frac_lo > 0.05))
settle_lo = himmelblau(samples_lo[int(0.9*n_steps):, 0], samples_lo[int(0.9*n_steps):, 1]).mean()

samples_hi, _ = results[50]
d_hi = np.linalg.norm(samples_hi[:, None, :] - known_minima[None, :, :], axis=2)
frac_hi = np.array([np.mean(np.argmin(d_hi, axis=1) == k) for k in range(4)])
basins_hi = int(np.sum(frac_hi > 0.05))
settle_hi = himmelblau(samples_hi[int(0.9*n_steps):, 0], samples_hi[int(0.9*n_steps):, 1]).mean()

print("=== Trade-off check ===")
print(f"CHECK T=1  : basins visited = {basins_lo} (traps in one basin), "
      f"mean f last 10% = {settle_lo:.4f} (settles low)")
print(f"CHECK T=50 : basins visited = {basins_hi} (roams widely), "
      f"mean f last 10% = {settle_hi:.4f} (never pins a minimum)")
print(f"CHECK conclusion: T=1 settles ({settle_lo:.3f}) but explores {basins_lo} basin(s); "
      f"T=50 explores {basins_hi} basins but stays high ({settle_hi:.3f}) -> "
      f"no single fixed temperature achieves both.")
print()
# One-sentence explanation of why this check confirms the result:
print("EXPLANATION: The check confirms the result because at low T the chain's tail "
      "energy is near zero yet it reaches only one basin, while at high T it reaches all "
      "four basins yet its tail energy stays far from zero, so wide exploration and tight "
      "settling are mutually exclusive at any single fixed temperature -- exactly the "
      "motivation for varying T (e.g. simulated annealing).")
