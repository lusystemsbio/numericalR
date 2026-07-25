import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Himmelblau's function ----
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# ---- setup ----
np.random.seed(1)                       # seed 1
n_steps = int(1e4)                       # 1e4 steps
mix_frac = 0.9                           # fraction of moves that are position moves
step = 1.0                               # step size for position proposals

ladder = np.arange(5, 85, 5, dtype=float)          # 16-rung ladder: 5,10,...,80
weights = np.linspace(5.0, 1.0, ladder.size)       # weights decrease linearly 5 -> 1
n_rungs = ladder.size

# the four known minima of Himmelblau (all with f = 0)
minima = np.array([[ 3.0,        2.0      ],
                   [-2.805118,   3.131312 ],
                   [-3.779310,  -3.283186 ],
                   [ 3.584428,  -1.848126 ]])

# ---- state ----
x, y = 0.0, 0.0                          # start at (0, 0)
k = 0                                    # current temperature-ladder index (start hottest? use lowest rung index)
fx = f(x, y)

# storage
xs = np.empty(n_steps); ys = np.empty(n_steps)
Ts = np.empty(n_steps)

# ---- simulated tempering: one chain, temperature random-walks the ladder ----
for i in range(n_steps):
    T = ladder[k]
    if np.random.rand() < mix_frac:
        # ----- position move: Metropolis at current temperature T -----
        xp = x + step * np.random.randn()
        yp = y + step * np.random.randn()
        fp = f(xp, yp)
        # accept with min(1, exp(-(f'-f)/T))
        if np.random.rand() < np.exp(-(fp - fx) / T):
            x, y, fx = xp, yp, fp
    else:
        # ----- temperature move: random walk to a neighbouring rung -----
        kp = k + (1 if np.random.rand() < 0.5 else -1)
        if 0 <= kp < n_rungs:            # reject proposals off the ladder ends
            Tp = ladder[kp]
            # acceptance a = min(1, (c'/c) * exp(-f*(1/T' - 1/T)))
            a = (weights[kp] / weights[k]) * np.exp(-fx * (1.0 / Tp - 1.0 / T))
            if np.random.rand() < a:
                k = kp
    # record state after this step
    xs[i], ys[i], Ts[i] = x, y, ladder[k]

# ---- check: which basins were visited (assign each point to nearest minimum) ----
pts = np.column_stack([xs, ys])
d2 = ((pts[:, None, :] - minima[None, :, :])**2).sum(axis=2)  # squared dist to each min
nearest = d2.argmin(axis=1)
# a point is "settled in a basin" if it is genuinely close to that minimum (cold behaviour)
close = np.sqrt(d2.min(axis=1)) < 0.75
basins_visited = np.unique(nearest[close])

# ---- printed numerical results ----
print(f"Ladder rungs: {ladder.tolist()}")
print(f"Ladder weights: {np.round(weights,4).tolist()}")
print(f"Total steps: {n_steps}")
print(f"Start point: (0.0, 0.0)")
print(f"Final point: ({x:.6f}, {y:.6f})")
print(f"Final f-value: {fx:.6e}")
print(f"Minimum f-value reached along path: {f(xs, ys).min():.6e}")
print(f"Temperature min / max visited: {Ts.min():.1f} / {Ts.max():.1f}")
print(f"Mean temperature over run: {Ts.mean():.4f}")
for j in range(4):
    print(f"Basin {j} at ({minima[j,0]:.4f}, {minima[j,1]:.4f}) settled-visits: {np.sum(close & (nearest==j))}")
print(f"Number of distinct basins settled into: {basins_visited.size}")
print(f"All four basins visited: {basins_visited.size == 4}")

# ---- figures: path over landscape + temperature over steps ----
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)

fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# left: single chain's path over the Himmelblau landscape
cs = ax[0].contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
fig.colorbar(cs, ax=ax[0], label="log(1 + f)")
ax[0].plot(xs, ys, color="white", lw=0.4, alpha=0.5, label="chain path")
ax[0].scatter(minima[:, 0], minima[:, 1], c="red", s=80, marker="*",
              edgecolors="k", zorder=5, label="four minima")
ax[0].scatter([0], [0], c="cyan", s=60, marker="o", edgecolors="k", zorder=5, label="start")
ax[0].set_title("Single chain path over Himmelblau landscape")
ax[0].set_xlabel("x"); ax[0].set_ylabel("y"); ax[0].legend(loc="upper right")

# right: temperature over steps (the random walk on the ladder)
ax[1].plot(np.arange(n_steps), Ts, lw=0.6, color="darkred")
ax[1].set_title("Temperature over steps (wandering the ladder)")
ax[1].set_xlabel("step"); ax[1].set_ylabel("temperature T")
ax[1].set_yticks(ladder)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.5.1_s4.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: Seeing settled visits (points within 0.75 of a minimum) accumulate in "
      "all four basins from a single unbroken trajectory confirms simulated tempering worked, "
      "because reaching every basin requires the chain to heat up and cross barriers between "
      "cold-phase descents into distinct minima.")
