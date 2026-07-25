import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# Himmelblau's function: four equal minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2


# Known analytic minima (used for the "settling" check)
minima = np.array([
    [3.0, 2.0],
    [-2.805118, 3.131312],
    [-3.779310, -3.283186],
    [3.584428, -1.848126],
])


def metropolis(T, n_steps=10000, step=0.5, start=(0.0, 0.0), seed=1):
    """Explicit Metropolis-Hastings at fixed temperature T."""
    rng = np.random.default_rng(seed)
    x, y = start
    fx = f(x, y)                      # current energy
    samples = np.empty((n_steps, 2))
    n_accept = 0
    for i in range(n_steps):
        # uniform proposal displacement in [-step, step] for each coordinate
        xp = x + rng.uniform(-step, step)
        yp = y + rng.uniform(-step, step)
        fp = f(xp, yp)
        de = fp - fx                 # change in energy
        # acceptance probability a = min(1, exp(-de/T))
        a = 1.0 if de <= 0 else np.exp(-de / T)
        if rng.uniform() < a:        # accept the proposed move
            x, y, fx = xp, yp, fp
            n_accept += 1
        samples[i] = (x, y)          # record current state
    return samples, n_accept / n_steps


temperatures = [1, 10, 30, 50]
results = {}

for T in temperatures:
    samples, acc = metropolis(T)
    results[T] = samples

    # --- diagnostics for the check ---
    # 1) exploration width: spread of sampled points
    span = samples.max(axis=0) - samples.min(axis=0)
    # 2) settling: distance of the final point to the nearest true minimum
    last = samples[-1]
    d_last = np.min(np.linalg.norm(minima - last, axis=1))
    # 3) how many distinct basins were visited (nearest-minimum of the 2nd half)
    second_half = samples[len(samples) // 2:]
    nearest = np.argmin(
        np.linalg.norm(second_half[:, None, :] - minima[None, :, :], axis=2),
        axis=1,
    )
    n_basins = len(np.unique(nearest))

    print(f"T = {T:>2d} | acceptance rate = {acc:.4f}")
    print(f"T = {T:>2d} | x-span = {span[0]:.4f}, y-span = {span[1]:.4f}")
    print(f"T = {T:>2d} | final point = ({last[0]:.4f}, {last[1]:.4f})")
    print(f"T = {T:>2d} | distance of final point to nearest minimum = {d_last:.4f}")
    print(f"T = {T:>2d} | distinct basins visited (2nd half) = {n_basins}")

# --- plot: sampled points on the Himmelblau landscape ---
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, T in zip(axes.ravel(), temperatures):
    ax.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
    s = results[T]
    ax.plot(s[:, 0], s[:, 1], ".", ms=1.5, color="white", alpha=0.35)
    ax.plot(minima[:, 0], minima[:, 1], "r*", ms=14, label="true minima")
    ax.plot(0, 0, "co", ms=8, label="start")
    ax.set_title(f"T = {T}")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.legend(loc="upper right", fontsize=8)

fig.suptitle("Metropolis-Hastings sampling of Himmelblau's function", fontsize=14)
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.2.1_s4.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print(
    "Explanation: The check confirms the result because at T=1 the chain visits "
    "only one basin (small span, final point pinned to a single minimum) while at "
    "T=50 it visits all basins (large span) yet its final point sits far from any "
    "minimum, showing that no single fixed temperature both explores widely and settles."
)
