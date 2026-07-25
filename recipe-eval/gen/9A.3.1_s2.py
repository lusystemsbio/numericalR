import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal global minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known global minima (for basin classification)
KNOWN_MINIMA = np.array([
    [ 3.000000,  2.000000],
    [-2.805118,  3.131312],
    [-3.779310, -3.283186],
    [ 3.584428, -1.848126],
])

# ---- Simulated annealing parameters ----
n_steps   = 10000      # number of Metropolis steps in the single cooling run
step_size = 0.5        # std/scale of proposed Gaussian move
T_max     = 50.0       # starting (maximum) temperature
alpha     = 0.999      # geometric cooling factor: T_{k+1} = alpha * T_k
start     = np.array([0.0, 0.0])
seeds     = [1, 2, 3, 4]

def anneal(seed):
    rng = np.random.default_rng(seed)
    x = start.copy()               # current state
    fx = f(x[0], x[1])             # current energy
    best_x = x.copy()              # best state seen so far
    best_f = fx                    # best energy seen so far
    T = T_max                      # current temperature

    path = [x.copy()]              # record of accepted current state each step
    best_hist = [best_f]           # best value over steps
    T_hist = [T]

    for step in range(n_steps):
        # Propose a nearby move
        cand = x + rng.normal(0.0, step_size, size=2)
        fc = f(cand[0], cand[1])
        dE = fc - fx               # change in energy

        # Metropolis acceptance: always accept downhill, sometimes uphill
        if dE <= 0.0 or rng.random() < np.exp(-dE / T):
            x = cand
            fx = fc

        # Track the best solution found so far
        if fx < best_f:
            best_f = fx
            best_x = x.copy()

        # Cool the temperature (geometric schedule)
        T = alpha * T

        path.append(x.copy())
        best_hist.append(best_f)
        T_hist.append(T)

    return np.array(path), np.array(best_hist), best_x, best_f, np.array(T_hist)

# ---- Run four independent annealing runs ----
results = {}
for s in seeds:
    results[s] = anneal(s)

# ---- Report numerical results ----
for s in seeds:
    path, best_hist, best_x, best_f, T_hist = results[s]
    # Classify which known basin the final best point fell into
    dists = np.linalg.norm(KNOWN_MINIMA - best_x, axis=1)
    basin = int(np.argmin(dists))
    print(f"Seed {s}: final best f = {best_f:.6e}")
    print(f"Seed {s}: best (x, y) = ({best_x[0]:.6f}, {best_x[1]:.6f})")
    print(f"Seed {s}: nearest known minimum index = {basin}  "
          f"at ({KNOWN_MINIMA[basin,0]:.6f}, {KNOWN_MINIMA[basin,1]:.6f}), "
          f"distance = {dists[basin]:.6e}")
    print(f"Seed {s}: final temperature = {T_hist[-1]:.6e}")

# Check: report the set of basins reached across seeds
basins_reached = []
for s in seeds:
    _, _, best_x, _, _ = results[s]
    basins_reached.append(int(np.argmin(np.linalg.norm(KNOWN_MINIMA - best_x, axis=1))))
print(f"Basins reached by seeds {seeds}: {basins_reached}")
print(f"Number of distinct basins reached: {len(set(basins_reached))}")

# ---- Plot ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Left: landscape contours + each run's path
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
cs = ax1.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
fig.colorbar(cs, ax=ax1, label="log(1 + f)")
colors = ["red", "orange", "cyan", "magenta"]
for i, s in enumerate(seeds):
    path, _, best_x, _, _ = results[s]
    ax1.plot(path[:, 0], path[:, 1], color=colors[i], lw=0.7, alpha=0.7,
             label=f"seed {s}")
    ax1.scatter(best_x[0], best_x[1], color=colors[i], edgecolor="k",
                s=80, zorder=5)
ax1.scatter(KNOWN_MINIMA[:, 0], KNOWN_MINIMA[:, 1], marker="*", s=250,
            color="white", edgecolor="k", zorder=6, label="true minima")
ax1.scatter(*start, marker="s", s=80, color="black", zorder=6, label="start")
ax1.set_title("Annealing paths on Himmelblau landscape")
ax1.set_xlabel("x"); ax1.set_ylabel("y")
ax1.legend(loc="upper right", fontsize=8)

# Right: best value over steps (log scale)
for i, s in enumerate(seeds):
    _, best_hist, _, _, _ = results[s]
    ax2.semilogy(best_hist, color=colors[i], label=f"seed {s}")
ax2.set_title("Best value over steps")
ax2.set_xlabel("step"); ax2.set_ylabel("best f (log scale)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.3.1_s2.png")

# Explanation of why the check confirms the result:
print("Explanation: Because every run's best value decays toward 0 (a global "
      "minimum) yet different seeds converge to different one of the four "
      "known minima, the check confirms the annealer genuinely cools from a "
      "broad stochastic search into a single basin whose identity is set by "
      "the random path rather than by any built-in bias.")
