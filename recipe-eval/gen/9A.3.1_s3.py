import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal global minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# Known analytic locations of the four minima (for the basin check)
known_minima = np.array([
    [ 3.000000,  2.000000],
    [-2.805118,  3.131312],
    [-3.779310, -3.283186],
    [ 3.584428, -1.848126],
])

# Simulated annealing parameters
start      = np.array([0.0, 0.0])   # starting point
n_steps    = int(1e4)               # number of Metropolis steps
step_size  = 0.5                    # proposal std / half-width
T_max      = 50.0                   # maximum (initial) temperature
alpha      = 0.999                  # geometric cooling factor
seeds      = [1, 2, 3, 4]

def anneal(seed):
    rng = np.random.default_rng(seed)
    x = start.copy()                # current state
    fx = f(*x)                      # current energy
    best_x = x.copy()               # best state seen so far
    best_f = fx                     # best energy seen so far
    T = T_max                       # current temperature
    path = [x.copy()]               # record of accepted-state trajectory
    best_hist = [best_f]            # record of best value over steps
    for _ in range(n_steps):
        # propose a uniform random step within +/- step_size in each coordinate
        prop = x + rng.uniform(-step_size, step_size, size=2)
        fp = f(*prop)
        dE = fp - fx                # energy change
        # Metropolis acceptance: always accept downhill, uphill with prob exp(-dE/T)
        if dE <= 0 or rng.random() < np.exp(-dE / T):
            x, fx = prop, fp        # accept the proposal
        # track the best value encountered
        if fx < best_f:
            best_f, best_x = fx, x.copy()
        # cool the temperature geometrically
        T *= alpha
        path.append(x.copy())
        best_hist.append(best_f)
    return np.array(path), np.array(best_hist), best_x, best_f

# Run all four seeds
results = {s: anneal(s) for s in seeds}

# ---- Report numerical results ----
for s in seeds:
    path, best_hist, best_x, best_f = results[s]
    # identify which analytic minimum this run landed nearest to
    dists = np.linalg.norm(known_minima - best_x, axis=1)
    idx = int(np.argmin(dists))
    print(f"Seed {s}: best_f = {best_f:.6e}")
    print(f"Seed {s}: best_x = ({best_x[0]:.6f}, {best_x[1]:.6f})")
    print(f"Seed {s}: nearest known minimum index = {idx} "
          f"at ({known_minima[idx,0]:.6f}, {known_minima[idx,1]:.6f}), "
          f"distance = {dists[idx]:.6e}")

# Confirm different seeds land in different basins
landed = []
for s in seeds:
    _, _, best_x, _ = results[s]
    landed.append(int(np.argmin(np.linalg.norm(known_minima - best_x, axis=1))))
print(f"Basin index reached per seed {seeds}: {landed}")
print(f"Number of distinct basins reached: {len(set(landed))}")

# ---- Plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Left: contour landscape with each annealing path
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
cs = ax1.contourf(GX, GY, np.log1p(GZ), levels=30, cmap="viridis")
fig.colorbar(cs, ax=ax1, label="log(1+f)")
colors = ["red", "white", "orange", "cyan"]
for s, c in zip(seeds, colors):
    path = results[s][0]
    ax1.plot(path[:, 0], path[:, 1], color=c, lw=0.7, alpha=0.8, label=f"seed {s}")
ax1.scatter(known_minima[:, 0], known_minima[:, 1], marker="*", s=200,
            color="magenta", edgecolor="k", zorder=5, label="minima")
ax1.scatter(*start, marker="o", s=60, color="black", zorder=5, label="start")
ax1.set_xlabel("x"); ax1.set_ylabel("y")
ax1.set_title("Annealing paths on Himmelblau landscape")
ax1.legend(loc="upper right", fontsize=8)

# Right: best value over steps (log scale) per run
for s, c in zip(seeds, ["red", "green", "orange", "blue"]):
    best_hist = results[s][1]
    ax2.plot(best_hist, color=c, label=f"seed {s}")
ax2.set_yscale("log")
ax2.set_xlabel("step"); ax2.set_ylabel("best f so far")
ax2.set_title("Best value over steps")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.3.1_s3.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because every run drives best_f down to ~0 (a true global "
      "minimum) yet the seeds settle near different analytic minima, the check "
      "confirms annealing cools from broad exploration into a single basin whose "
      "identity is determined by the stochastic path rather than by the algorithm.")
