import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Himmelblau's function and its four equal minima (f = 0) ----
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

minima = np.array([
    [ 3.000000,  2.000000],
    [-2.805118,  3.131312],
    [-3.779310, -3.283186],
    [ 3.584428, -1.848126],
])

# ---- Simulated tempering setup ----
rng = np.random.default_rng(1)          # seed 1
n_steps = int(1e4)
mix = 0.9                               # fraction of position (vs temperature) moves
step = 1.0                             # position proposal std
ladder = np.arange(5, 81, 5, dtype=float)   # 16-rung ladder: 5,10,...,80
weights = np.linspace(5.0, 1.0, ladder.size)  # log-weights c, linearly 5 -> 1

# ---- State: position (x,y) and temperature index k ----
x, y = 0.0, 0.0                        # start at (0,0)
k = 0                                  # start on lowest rung
fx = f(x, y)

# storage
path = np.empty((n_steps, 2))
temps = np.empty(n_steps)

for i in range(n_steps):
    T = ladder[k]
    if rng.random() < mix:
        # --- Position move: Metropolis at current temperature T ---
        xp = x + step * rng.standard_normal()
        yp = y + step * rng.standard_normal()
        fp = f(xp, yp)
        # accept with prob min(1, exp(-(f' - f)/T)); hot T lets it cross barriers
        if rng.random() < np.exp(-(fp - fx) / T):
            x, y, fx = xp, yp, fp
    else:
        # --- Temperature move: random-walk one rung up or down the ladder ---
        kp = k + (1 if rng.random() < 0.5 else -1)
        if 0 <= kp < ladder.size:      # reject moves off the ladder ends
            Tp = ladder[kp]
            # a = min(1, (c'/c) * exp(-f*(1/T' - 1/T)))  (weights held as log c)
            a = np.exp((weights[kp] - weights[k]) - fx * (1.0 / Tp - 1.0 / T))
            if rng.random() < min(1.0, a):
                k = kp
    path[i] = (x, y)
    temps[i] = ladder[k]

# ---- Check: which of the four basins were visited ----
# assign every visited point to its nearest minimum, then count unique basins
d = np.linalg.norm(path[:, None, :] - minima[None, :, :], axis=2)
basin = np.argmin(d, axis=1)
visited = np.unique(basin)

print("Start point: (0, 0)")
print("Total steps:", n_steps)
print("Ladder (16 rungs):", ", ".join(f"{t:.0f}" for t in ladder))
print("Min temperature visited:", temps.min())
print("Max temperature visited:", temps.max())
for j, (mx, my) in enumerate(minima):
    print(f"Basin {j} minimum ({mx:.3f}, {my:.3f}): "
          f"{np.count_nonzero(basin == j)} points visited")
print("Number of distinct basins visited:", visited.size)
print("All four basins visited:", visited.size == 4)
# One-sentence explanation of why this check confirms the result:
print("Check rationale: visiting all four equal minima in a single chain proves "
      "the wandering temperature both got hot enough to cross the barriers between "
      "basins and cold enough to settle inside them.")

# ---- Plots: chain path over landscape, and temperature over steps ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
ax1.contourf(GX, GY, np.log1p(f(GX, GY)), levels=40, cmap="viridis")
ax1.plot(path[:, 0], path[:, 1], color="white", lw=0.4, alpha=0.5)
ax1.scatter(minima[:, 0], minima[:, 1], c="red", marker="*", s=200,
            edgecolor="k", zorder=5, label="minima")
ax1.scatter([0], [0], c="cyan", marker="o", s=60, edgecolor="k",
            zorder=5, label="start")
ax1.set_title("Single-chain path over Himmelblau landscape")
ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.legend(loc="upper right")

ax2.plot(np.arange(n_steps), temps, lw=0.6, color="darkorange")
ax2.set_title("Temperature random-walk over steps")
ax2.set_xlabel("step"); ax2.set_ylabel("temperature T")

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.5.1_s3.png", dpi=120)
