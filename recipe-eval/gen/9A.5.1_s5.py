import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Himmelblau's function: four equal minima at f = 0 ---
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known minima of Himmelblau's function
true_minima = np.array([
    [ 3.0,        2.0],
    [-2.805118,   3.131312],
    [-3.779310,  -3.283186],
    [ 3.584428,  -1.848126],
])

# --- Setup ---
rng = np.random.default_rng(1)          # seed 1
n_steps = int(1e4)                       # 1e4 steps
mix_frac = 0.9                           # mixing fraction: prob of a spatial (vs temperature) move
step_size = 1.0                          # spatial proposal step size

# 16-rung temperature ladder T = 5, 10, ..., 80
ladder = np.arange(5, 85, 5, dtype=float)
n_rungs = len(ladder)

# Weights (pseudo-counts c) decreasing linearly from 5 to 1 across the ladder
weights = np.linspace(5.0, 1.0, n_rungs)

# --- Chain state ---
x, y = 0.0, 0.0                          # start at (0, 0)
k = 0                                    # start on the lowest rung (T = 5)

# Storage
xs = np.empty(n_steps + 1)
ys = np.empty(n_steps + 1)
Ts = np.empty(n_steps + 1)
xs[0], ys[0], Ts[0] = x, y, ladder[k]

fx = f(x, y)                             # cache current function value

# --- Simulated tempering: single chain, temperature random-walks the ladder ---
for i in range(1, n_steps + 1):
    T = ladder[k]
    if rng.random() < mix_frac:
        # --- Spatial Metropolis move at current temperature T ---
        xp = x + rng.normal(0.0, step_size)
        yp = y + rng.normal(0.0, step_size)
        fp = f(xp, yp)
        # accept with prob min(1, exp(-(f' - f)/T))
        if rng.random() < np.exp(-(fp - fx) / T):
            x, y, fx = xp, yp, fp
    else:
        # --- Temperature move: random-walk to a neighboring rung ---
        kp = k + (1 if rng.random() < 0.5 else -1)
        if 0 <= kp < n_rungs:
            Tp = ladder[kp]
            cp, c = weights[kp], weights[k]     # target vs current rung weight
            # accept with prob a = min(1, (c'/c) * exp(-f*(1/T' - 1/T)))
            a = (cp / c) * np.exp(-fx * (1.0 / Tp - 1.0 / T))
            if rng.random() < a:
                k = kp
    xs[i], ys[i], Ts[i] = x, y, ladder[k]

# --- Check: assign each visited point to its nearest true minimum basin ---
# Consider a point "settled" in a basin if it is cold (low T) and close to a minimum.
settled = Ts <= 10.0                                  # cold portion of the trajectory
pts = np.column_stack([xs, ys])
# nearest-minimum index for each point
dists = np.linalg.norm(pts[:, None, :] - true_minima[None, :, :], axis=2)
nearest = np.argmin(dists, axis=1)
near_dist = np.min(dists, axis=1)
in_basin = settled & (near_dist < 0.5)                # cold AND near a minimum
visited_basins = sorted(set(nearest[in_basin].tolist()))

# --- Numerical results ---
print(f"Number of steps: {n_steps}")
print(f"Number of temperature rungs: {n_rungs}")
print(f"Temperature ladder: {ladder.tolist()}")
print(f"Rung weights: {weights.tolist()}")
print(f"Min temperature visited: {Ts.min()}")
print(f"Max temperature visited: {Ts.max()}")
print(f"Mean temperature: {Ts.mean():.4f}")
print(f"Final function value f(x,y): {fx:.6f}")
print(f"Best (lowest) function value along path: {f(xs, ys).min():.6f}")
for j, (mx, my) in enumerate(true_minima):
    cnt = int(np.sum(in_basin & (nearest == j)))
    print(f"Basin {j} at ({mx:.4f}, {my:.4f}): cold visits = {cnt}")
print(f"Number of distinct basins visited (cold & near-minimum): {len(visited_basins)}")
print(f"All four basins visited: {len(visited_basins) == 4}")

# --- Plots ---
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# (1) Chain path over the Himmelblau landscape
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
cs = axes[0].contourf(GX, GY, np.log1p(GZ), levels=30, cmap="viridis")
fig.colorbar(cs, ax=axes[0], label="log(1 + f)")
axes[0].plot(xs, ys, color="white", lw=0.4, alpha=0.6)
axes[0].scatter(true_minima[:, 0], true_minima[:, 1], c="red", s=120,
                marker="*", edgecolor="black", zorder=5, label="true minima")
axes[0].scatter([0], [0], c="cyan", s=60, marker="o", zorder=5, label="start (0,0)")
axes[0].set_title("Single chain path over Himmelblau landscape")
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
axes[0].legend(loc="upper right")

# (2) Temperature over steps
axes[1].plot(np.arange(n_steps + 1), Ts, color="darkorange", lw=0.6)
axes[1].set_title("Temperature random-walk over steps")
axes[1].set_xlabel("step"); axes[1].set_ylabel("temperature T")
axes[1].set_yticks(ladder)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.5.1_s5.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: Finding cold visits clustered near all four distinct minima in a "
      "single unbroken chain confirms that letting the temperature wander lets one "
      "trajectory hop barriers while hot and then relax into each basin when cold.")
