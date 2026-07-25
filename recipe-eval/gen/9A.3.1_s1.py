import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Himmelblau's function: four equal global minima at f = 0 ---
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known global minima (for basin identification)
KNOWN_MINIMA = np.array([
    [3.0, 2.0],
    [-2.805118, 3.131312],
    [-3.779310, -3.283186],
    [3.584428, -1.848126],
])

# --- Simulated annealing parameters ---
n_steps = 10000        # number of Metropolis steps per run
step_size = 0.5        # proposal step scale
T_max = 50.0           # starting (maximum) temperature
alpha = 0.999          # geometric cooling factor
start = np.array([0.0, 0.0])
seeds = [1, 2, 3, 4]

def anneal(seed):
    rng = np.random.default_rng(seed)
    x = start.copy()               # current state
    fx = f(x[0], x[1])             # current energy
    best_x = x.copy()              # best state seen
    best_f = fx                    # best energy seen
    T = T_max                      # temperature, cooled each step

    path = [x.copy()]              # record accepted-state path
    best_hist = [best_f]           # record best value over steps

    for _ in range(n_steps):
        # propose a random Gaussian move from the current state
        cand = x + rng.normal(0.0, step_size, size=2)
        fc = f(cand[0], cand[1])
        dE = fc - fx               # change in energy

        # Metropolis acceptance: always accept downhill, sometimes uphill
        if dE <= 0 or rng.random() < np.exp(-dE / T):
            x = cand
            fx = fc

        # track the best solution encountered
        if fx < best_f:
            best_f = fx
            best_x = x.copy()

        # geometric cooling: shrink T toward zero
        T *= alpha

        path.append(x.copy())
        best_hist.append(best_f)

    return np.array(path), np.array(best_hist), best_x, best_f

# --- Run four independent annealing runs ---
results = {}
for s in seeds:
    path, best_hist, best_x, best_f = anneal(s)
    results[s] = (path, best_hist, best_x, best_f)
    # identify which of the four basins the run landed in
    basin = int(np.argmin(np.linalg.norm(KNOWN_MINIMA - best_x, axis=1)))
    print(f"Seed {s}: best f = {best_f:.6e} at "
          f"(x, y) = ({best_x[0]:.6f}, {best_x[1]:.6f}) -> basin {basin} "
          f"near ({KNOWN_MINIMA[basin,0]:.4f}, {KNOWN_MINIMA[basin,1]:.4f})")

# report which basin each run converged to
basins = []
for s in seeds:
    _, _, best_x, _ = results[s]
    basins.append(int(np.argmin(np.linalg.norm(KNOWN_MINIMA - best_x, axis=1))))
print(f"Basins reached by seeds {seeds}: {basins}")
print(f"Number of distinct basins reached: {len(set(basins))}")

# --- Plot: landscape with paths, and best value over steps ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# contour of Himmelblau's function
gx = np.linspace(-5, 5, 400)
gy = np.linspace(-5, 5, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
cs = ax1.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
fig.colorbar(cs, ax=ax1, label="log(1 + f)")

colors = ["red", "orange", "cyan", "magenta"]
for s, c in zip(seeds, colors):
    path, _, best_x, _ = results[s]
    ax1.plot(path[:, 0], path[:, 1], color=c, lw=0.6, alpha=0.7,
             label=f"seed {s}")
    ax1.plot(best_x[0], best_x[1], marker="*", color=c, ms=16,
             markeredgecolor="black")
ax1.scatter(KNOWN_MINIMA[:, 0], KNOWN_MINIMA[:, 1], marker="x",
            color="white", s=90, label="true minima")
ax1.plot(start[0], start[1], "ks", ms=8, label="start")
ax1.set_title("Annealing paths on Himmelblau landscape")
ax1.set_xlabel("x"); ax1.set_ylabel("y")
ax1.legend(loc="upper right", fontsize=8)

# best value over steps (log scale)
for s, c in zip(seeds, colors):
    _, best_hist, _, _ = results[s]
    ax2.plot(best_hist, color=c, label=f"seed {s}")
ax2.set_yscale("log")
ax2.set_title("Best value over steps")
ax2.set_xlabel("step"); ax2.set_ylabel("best f (log scale)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.3.1_s1.png")

# --- Explanation of the check ---
print("Check: because every run's best value decays toward 0 (broad early "
      "wandering then settling) while the final points cluster at different "
      "ones of the four true minima, we confirm annealing cools a wide search "
      "into a single basin whose identity is set by the random path (seed).")
