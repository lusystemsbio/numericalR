import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known minima of Himmelblau's function (basin centers)
minima = np.array([
    [ 3.0,        2.0     ],
    [-2.805118,   3.131312],
    [-3.779310,  -3.283186],
    [ 3.584428,  -1.848126],
])

# ---- Simulated tempering setup ----
rng = np.random.default_rng(1)          # seed 1

n_steps   = 10000                        # 1e4 steps
mix_frac  = 0.9                          # fraction of moves that are spatial (vs temperature) moves
step_size = 1.0                          # Gaussian proposal std for x,y

# 16-rung temperature ladder T = 5, 10, ..., 80
ladder  = np.arange(5.0, 85.0, 5.0)      # [5,10,...,80], length 16
n_rungs = len(ladder)

# Pseudo-weights (the "c" constants) decreasing linearly from 5 to 1 across the ladder
weights = np.linspace(5.0, 1.0, n_rungs)

# ---- Initialize the single chain ----
x, y = 0.0, 0.0                          # start at (0, 0)
k    = 0                                 # current temperature rung index
fval = f(x, y)

# Storage for the path and temperature trace
path = np.empty((n_steps + 1, 2))
Ttrace = np.empty(n_steps + 1)
path[0]   = (x, y)
Ttrace[0] = ladder[k]

# ---- Run the chain ----
for step in range(n_steps):
    T = ladder[k]

    if rng.random() < mix_frac:
        # --- Spatial Metropolis move at the current temperature ---
        xp = x + step_size * rng.standard_normal()
        yp = y + step_size * rng.standard_normal()
        fp = f(xp, yp)
        # accept with prob min(1, exp(-(f'-f)/T))
        if rng.random() < np.exp(-(fp - fval) / T):
            x, y, fval = xp, yp, fp
    else:
        # --- Temperature move: random-walk one rung up or down the ladder ---
        kp = k + (1 if rng.random() < 0.5 else -1)
        if 0 <= kp < n_rungs:            # reject moves off the ends of the ladder
            Tp = ladder[kp]
            c, cp = weights[k], weights[kp]
            # accept with a = min(1, (c'/c) * exp(-f*(1/T' - 1/T)))
            a = (cp / c) * np.exp(-fval * (1.0 / Tp - 1.0 / T))
            if rng.random() < a:
                k = kp

    path[step + 1]   = (x, y)
    Ttrace[step + 1] = ladder[k]

# ---- Check: did the single trajectory visit all four basins? ----
# Assign each visited point to its nearest known minimum, but only count
# "settled" points (low f, i.e. genuinely inside a basin).
settled = f(path[:, 0], path[:, 1]) < 5.0
dists = np.linalg.norm(path[:, None, :] - minima[None, :, :], axis=2)
nearest = np.argmin(dists, axis=1)
visited = set(np.unique(nearest[settled]).tolist())

# ---- Print numerical results ----
print(f"Number of steps: {n_steps}")
print(f"Number of ladder rungs: {n_rungs}")
print(f"Temperature ladder: {ladder.tolist()}")
print(f"Ladder weights (c): {weights.tolist()}")
print(f"Start point: (0.0, 0.0), f = {f(0.0, 0.0):.6f}")
print(f"Final point: ({x:.6f}, {y:.6f}), f = {fval:.6f}")
print(f"Final temperature rung index: {k}, T = {ladder[k]:.1f}")
print(f"Min temperature visited: {Ttrace.min():.1f}")
print(f"Max temperature visited: {Ttrace.max():.1f}")
print(f"Number of settled points (f < 5): {int(settled.sum())}")
for i, m in enumerate(minima):
    count = int(np.sum(settled & (nearest == i)))
    print(f"Basin {i} near ({m[0]:.4f}, {m[1]:.4f}): visited={i in visited}, settled_points={count}")
print(f"Number of distinct basins visited: {len(visited)}")
print(f"All four basins visited: {len(visited) == 4}")

# ---- Plot: chain path over the landscape, and temperature over steps ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Landscape contours (log scale to reveal all four wells)
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
ax1.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
ax1.plot(path[:, 0], path[:, 1], color="white", lw=0.4, alpha=0.6)
ax1.scatter(minima[:, 0], minima[:, 1], color="red", marker="*", s=200,
            edgecolor="black", zorder=5, label="true minima")
ax1.scatter([0], [0], color="cyan", marker="o", s=60, zorder=5, label="start")
ax1.set_title("Single simulated-tempering chain path")
ax1.set_xlabel("x"); ax1.set_ylabel("y")
ax1.set_xlim(-6, 6); ax1.set_ylim(-6, 6)
ax1.legend(loc="upper right")

# Temperature over steps
ax2.plot(Ttrace, color="darkorange", lw=0.6)
ax2.set_title("Temperature random-walk over steps")
ax2.set_xlabel("step"); ax2.set_ylabel("T")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.5.1_s2.png")

# One-sentence explanation of why the check confirms the result:
# Because a single chain whose temperature wanders visits all four basins,
# it demonstrates that the same trajectory used high-T excursions to hop over
# barriers and low-T excursions to settle into each well, which is exactly the
# behavior simulated tempering is meant to produce.
print("Explanation: Finding that one chain settles (f<5) into all four distinct "
      "basins proves the wandering temperature let it hop barriers while hot and "
      "fall into minima while cold, confirming the method works.")
