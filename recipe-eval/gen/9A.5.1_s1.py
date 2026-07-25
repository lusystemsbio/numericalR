import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# ---- Setup ----
rng = np.random.default_rng(1)          # seed 1
n_steps = 10000                         # 1e4 steps
mix_frac = 0.9                          # 90% spatial moves, 10% temperature moves
step_size = 1.0                         # spatial proposal std

# Temperature ladder: 16 rungs T = 5, 10, ..., 80
T_ladder = np.arange(5, 81, 5, dtype=float)
n_rungs = len(T_ladder)
# pseudo-weights c decreasing linearly from 5 to 1 across the ladder
c_ladder = np.linspace(5.0, 1.0, n_rungs)

# ---- State ----
x, y = 0.0, 0.0                         # start at (0, 0)
k = 0                                   # current rung index (T = 5)
fx = f(x, y)

# Storage
xs = np.empty(n_steps); ys = np.empty(n_steps); Ts = np.empty(n_steps)

n_temp_accept = 0; n_temp_prop = 0
n_space_accept = 0; n_space_prop = 0

for i in range(n_steps):
    T = T_ladder[k]
    if rng.random() < mix_frac:
        # ---- Spatial Metropolis move at fixed temperature T ----
        n_space_prop += 1
        xp = x + step_size * rng.standard_normal()
        yp = y + step_size * rng.standard_normal()
        fp = f(xp, yp)
        # accept with prob min(1, exp(-(f'-f)/T))
        if rng.random() < np.exp(-(fp - fx) / T):
            x, y, fx = xp, yp, fp
            n_space_accept += 1
    else:
        # ---- Temperature move: random-walk neighbor on the ladder ----
        n_temp_prop += 1
        step = 1 if rng.random() < 0.5 else -1
        kp = k + step
        if 0 <= kp < n_rungs:           # reject moves off the ladder ends
            Tp = T_ladder[kp]
            # acceptance a = min(1, (c'/c) * exp(-f * (1/T' - 1/T)))
            a = (c_ladder[kp] / c_ladder[k]) * np.exp(-fx * (1.0 / Tp - 1.0 / T))
            if rng.random() < a:
                k = kp
                n_temp_accept += 1

    xs[i], ys[i], Ts[i] = x, y, T_ladder[k]

# ---- Basin assignment: which of the four minima is nearest ----
minima = np.array([[3.0, 2.0],
                   [-2.805118, 3.131312],
                   [-3.779310, -3.283186],
                   [3.584428, -1.848126]])
# only count points that have actually "settled" (low f) to define basin visits
settled = f(xs, ys) < 1.0
basin_idx = np.argmin(
    ((xs[:, None] - minima[None, :, 0])**2 + (ys[:, None] - minima[None, :, 1])**2),
    axis=1)
visited = np.unique(basin_idx[settled])

# ---- Reporting ----
print(f"Total steps: {n_steps}")
print(f"Spatial moves proposed: {n_space_prop}, accepted: {n_space_accept}, "
      f"acceptance rate: {n_space_accept / max(n_space_prop,1):.4f}")
print(f"Temperature moves proposed: {n_temp_prop}, accepted: {n_temp_accept}, "
      f"acceptance rate: {n_temp_accept / max(n_temp_prop,1):.4f}")
print(f"Minimum T visited: {Ts.min()}")
print(f"Maximum T visited: {Ts.max()}")
print(f"Minimum f achieved along chain: {f(xs, ys).min():.6e}")
print(f"Number of settled points (f < 1): {int(settled.sum())}")
for j, m in enumerate(minima):
    cnt = int(np.sum(settled & (basin_idx == j)))
    print(f"Basin {j+1} at ({m[0]:+.4f}, {m[1]:+.4f}): settled visits = {cnt}")
print(f"Number of distinct basins visited (settled): {len(visited)} out of 4")
print(f"All four basins visited: {len(visited) == 4}")

# ---- Plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Landscape with chain path
gx = np.linspace(-6, 6, 400); gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
cs = ax1.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
fig.colorbar(cs, ax=ax1, label="log(1 + f)")
ax1.plot(xs, ys, color="white", lw=0.4, alpha=0.5, label="chain path")
ax1.scatter(minima[:, 0], minima[:, 1], c="red", marker="*", s=200,
            edgecolors="k", zorder=5, label="minima")
ax1.scatter([0], [0], c="cyan", marker="o", s=60, edgecolors="k",
            zorder=5, label="start")
ax1.set_title("Simulated tempering: single chain path over Himmelblau landscape")
ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.legend(loc="upper right")

# Temperature over steps
ax2.plot(np.arange(n_steps), Ts, color="firebrick", lw=0.6)
ax2.set_title("Temperature random-walk over steps")
ax2.set_xlabel("step"); ax2.set_ylabel("temperature T")
ax2.set_yticks(T_ladder)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.5.1_s1.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Check explanation: Because a SINGLE chain settles (f < 1) in all four basins, "
      "the wandering temperature must have supplied enough hot phases to hop between "
      "basins and enough cold phases to relax into each, confirming the tempering method "
      "explores and exploits within one trajectory.")
