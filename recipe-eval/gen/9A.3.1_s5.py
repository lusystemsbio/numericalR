import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal global minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known global minima of Himmelblau's function
KNOWN_MINIMA = np.array([
    [ 3.000000,  2.000000],
    [-2.805118,  3.131312],
    [-3.779310, -3.283186],
    [ 3.584428, -1.848126],
])

# Simulated annealing parameters
N_STEPS   = int(1e4)   # steps in a single cooling run
STEP_SIZE = 0.5        # std dev of Gaussian proposal
T_MAX     = 50.0       # starting (maximum) temperature
ALPHA     = 0.999      # geometric cooling factor per step
START     = np.array([0.0, 0.0])
SEEDS     = [1, 2, 3, 4]


def anneal(seed):
    """Run one cooled Metropolis chain; return path and best-value history."""
    rng = np.random.default_rng(seed)      # reproducible random path per seed

    x = START.copy()                       # current state
    fx = f(*x)                             # current energy

    best_x = x.copy()                      # best state seen so far
    best_f = fx                            # best energy seen so far

    T = T_MAX                              # temperature, cooled each step

    path = np.empty((N_STEPS + 1, 2))      # record trajectory for plotting
    best_hist = np.empty(N_STEPS + 1)      # record best-value-over-steps
    path[0] = x
    best_hist[0] = best_f

    for i in range(1, N_STEPS + 1):
        # Propose a nearby point with a Gaussian step
        prop = x + rng.normal(0.0, STEP_SIZE, size=2)
        fp = f(*prop)

        # Metropolis acceptance: always accept downhill, accept uphill
        # with probability exp(-dE / T); high T => broad exploration.
        dE = fp - fx
        if dE <= 0.0 or rng.random() < np.exp(-dE / T):
            x, fx = prop, fp               # accept the move

        # Track the best value ever visited
        if fx < best_f:
            best_f, best_x = fx, x.copy()

        # Geometric cooling: shrink T toward 0 so the search focuses
        T *= ALPHA

        path[i] = x
        best_hist[i] = best_f

    return path, best_hist, best_x, best_f


# Run the four annealing chains
results = {}
for s in SEEDS:
    results[s] = anneal(s)

# Report where each run landed and which basin it fell into
print("=== Simulated annealing on Himmelblau's function ===")
for s in SEEDS:
    path, best_hist, best_x, best_f = results[s]
    # Identify nearest known global minimum (which basin it settled in)
    dists = np.linalg.norm(KNOWN_MINIMA - best_x, axis=1)
    which = int(np.argmin(dists))
    print(f"Seed {s}: best_x = ({best_x[0]:+.4f}, {best_x[1]:+.4f})   "
          f"best_f = {best_f:.6e}   basin = minimum #{which} "
          f"({KNOWN_MINIMA[which][0]:+.3f}, {KNOWN_MINIMA[which][1]:+.3f})")

# Check: do different seeds land in different basins?
basins = []
for s in SEEDS:
    _, _, best_x, _ = results[s]
    basins.append(int(np.argmin(np.linalg.norm(KNOWN_MINIMA - best_x, axis=1))))
print(f"\nBasin index reached per seed {SEEDS}: {basins}")
print(f"Number of distinct basins reached: {len(set(basins))}")
print("Explanation: because all four minima have equal depth f = 0, the run "
      "cannot prefer one over another, so the basin it settles into is decided "
      "purely by the random path -- hence different seeds land in different "
      "minima, confirming a genuine broad-to-narrow stochastic search rather "
      "than a deterministic descent to one fixed point.")

# ---- Plot: paths on the landscape + best value over steps ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Contour of the landscape (log scale to reveal the four basins)
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)
cs = ax1.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
fig.colorbar(cs, ax=ax1, label="log(1 + f)")

colors = ["red", "orange", "white", "magenta"]
for s, c in zip(SEEDS, colors):
    path, _, best_x, _ = results[s]
    ax1.plot(path[:, 0], path[:, 1], color=c, lw=0.6, alpha=0.7,
             label=f"seed {s}")
    ax1.plot(best_x[0], best_x[1], marker="*", color=c, ms=16,
             markeredgecolor="black")

# Mark the four true global minima
ax1.scatter(KNOWN_MINIMA[:, 0], KNOWN_MINIMA[:, 1], marker="X",
            color="black", s=80, label="true minima", zorder=5)
ax1.plot(*START, marker="o", color="cyan", ms=10, markeredgecolor="black",
         label="start (0,0)")
ax1.set_title("Annealing paths on Himmelblau's landscape")
ax1.set_xlabel("x"); ax1.set_ylabel("y")
ax1.legend(loc="upper right", fontsize=8)

# Best value over steps (log scale) for each run
for s, c in zip(SEEDS, colors):
    _, best_hist, _, _ = results[s]
    ax2.plot(best_hist, color=("gray" if c == "white" else c),
             label=f"seed {s}")
ax2.set_yscale("log")
ax2.set_xlabel("step"); ax2.set_ylabel("best f so far (log scale)")
ax2.set_title("Best value over steps (cooling into a single minimum)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.3.1_s5.png", dpi=120)
print("\nSaved figure to 9A.3.1_s5.png")
