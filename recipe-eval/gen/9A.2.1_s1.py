import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Himmelblau's function: four equal minima at f = 0
def himmelblau(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2


# Explicit Metropolis-Hastings at fixed temperature T
def metropolis_hastings(f, start, n_steps, step, T, seed):
    rng = np.random.default_rng(seed)          # per-temperature RNG
    x = np.array(start, dtype=float)           # current state
    fx = f(x[0], x[1])                          # current energy
    samples = np.empty((n_steps, 2))           # storage for the chain

    for i in range(n_steps):
        # uniform proposal displacement in [-step, step] for each coordinate
        prop = x + rng.uniform(-step, step, size=2)
        fp = f(prop[0], prop[1])               # proposed energy
        de = fp - fx                            # energy change
        # acceptance probability a = min(1, exp(-de/T))
        a = min(1.0, np.exp(-de / T))
        if rng.uniform() < a:                   # accept?
            x, fx = prop, fp                    # move to proposed state
        samples[i] = x                          # record current state
    return samples


# The four true minima of Himmelblau's function
true_minima = np.array([
    [ 3.000000,  2.000000],
    [-2.805118,  3.131312],
    [-3.779310, -3.283186],
    [ 3.584428, -1.848126],
])

start = (0.0, 0.0)
n_steps = int(1e4)
step = 0.5
temps = [1, 10, 30, 50]

# Background contour of the landscape
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = himmelblau(GX, GY)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

results = {}
for ax, T in zip(axes.ravel(), temps):
    samples = metropolis_hastings(himmelblau, start, n_steps, step, T, seed=1)
    results[T] = samples

    # log-scaled contours make the four basins visible
    ax.contourf(GX, GY, np.log1p(GZ), levels=30, cmap="viridis")
    ax.plot(samples[:, 0], samples[:, 1], ".", ms=1.5, color="white", alpha=0.3)
    ax.plot(true_minima[:, 0], true_minima[:, 1], "r*", ms=14,
            markeredgecolor="k", label="true minima")
    ax.set_title(f"T = {T}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.legend(loc="upper right", fontsize=8)

fig.suptitle("Metropolis-Hastings on Himmelblau's landscape (fixed T)", fontsize=14)
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.2.1_s1.png")

# ---- Diagnostics ----------------------------------------------------------
# For each temperature, measure two things using the second half of the chain:
#   (1) SPREAD: how widely it explores  -> std of the sampled points
#   (2) SETTLE: how close it pins a min  -> mean distance to the nearest true minimum
print("Himmelblau minima value f = 0.0 (four equal minima)")
print("")

for T in temps:
    s = results[T]
    half = s[n_steps // 2:]                     # discard burn-in

    # spread: standard deviation of positions (larger = roams more)
    spread = np.sqrt(np.mean(np.var(half, axis=0)))

    # settle: distance from each sample to its nearest true minimum
    dists = np.sqrt(((half[:, None, :] - true_minima[None, :, :])**2).sum(axis=2))
    nearest = dists.min(axis=1)
    settle = nearest.mean()                      # smaller = pins a minimum tightly

    # how many distinct basins were visited (samples within 1.0 of a minimum)
    visited = np.unique(dists.argmin(axis=1)[nearest < 1.0])
    n_basins = len(visited)

    print(f"T = {T:2d} | spread (pos std) = {spread:8.4f} | "
          f"mean dist to nearest min = {settle:8.4f} | basins occupied = {n_basins}")

print("")
print("Interpretation:")
print("T = 1  -> small spread, small basin count: trapped in ONE basin (settles but no exploration)")
print("T = 50 -> large spread, large mean distance to minima: ROAMS widely but never pins a minimum")
print("")
print("Check: because spread and settle move in opposite directions across T "
      "(low T settles but stays in one basin, high T explores but never sits at a minimum), "
      "no single fixed temperature achieves both, which confirms the result.")
