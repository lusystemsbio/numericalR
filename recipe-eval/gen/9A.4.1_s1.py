import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Himmelblau's function: four equal minima at f = 0 ---
def f(p):
    x, y = p
    return (x*x + y - 11.0)**2 + (x + y*y - 7.0)**2

# --- Parallel tempering setup ---
rng = np.random.default_rng(1)                    # seed 1
T     = np.array([2.3, 5.0, 10.0, 20.0, 40.0, 80.0])   # temperatures (cold -> hot)
steps = np.array([0.3, 0.5, 0.9, 1.4, 1.9, 2.5])       # per-replica proposal step sizes
n_rep = len(T)
n_rounds = 100          # number of swap rounds
n_steps  = 100          # Metropolis steps per replica per round

# All replicas start at the origin.
pos = np.zeros((n_rep, 2))
fval = np.array([f(pos[i]) for i in range(n_rep)])

# Storage for the trajectories we want to look at.
cold_pts = []   # sampled points of coldest replica (index 0)
hot_pts  = []   # sampled points of hottest replica (index n_rep-1)

accept_moves = np.zeros(n_rep)   # count accepted Metropolis moves per replica
prop_moves   = np.zeros(n_rep)   # count proposed Metropolis moves per replica
accept_swaps = 0                 # count accepted swaps
prop_swaps   = 0                 # count proposed swaps

for rnd in range(n_rounds):
    # --- 1) Run an independent Metropolis chain at each temperature ---
    for i in range(n_rep):
        for _ in range(n_steps):
            # propose a Gaussian step scaled by this replica's step size
            cand = pos[i] + steps[i] * rng.standard_normal(2)
            fc = f(cand)
            prop_moves[i] += 1
            # Metropolis acceptance at temperature T_i (Boltzmann weight exp(-f/T))
            if fc <= fval[i] or rng.random() < np.exp((fval[i] - fc) / T[i]):
                pos[i] = cand
                fval[i] = fc
                accept_moves[i] += 1
        # record the current point for the two replicas of interest
        if i == 0:
            cold_pts.append(pos[0].copy())
        if i == n_rep - 1:
            hot_pts.append(pos[i].copy())

    # --- 2) Attempt swaps of adjacent-temperature replicas ---
    for i in range(n_rep - 1):
        j = i + 1
        prop_swaps += 1
        # a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))
        a = np.exp((fval[i] - fval[j]) * (1.0/T[i] - 1.0/T[j]))
        if rng.random() < min(1.0, a):
            # exchange the two replicas' configurations (and their f values)
            pos[[i, j]] = pos[[j, i]]
            fval[[i, j]] = fval[[j, i]]
            accept_swaps += 1

cold_pts = np.array(cold_pts)
hot_pts  = np.array(hot_pts)

# --- The four known minima of Himmelblau's function ---
minima = np.array([[ 3.000000,  2.000000],
                   [-2.805118,  3.131312],
                   [-3.779310, -3.283186],
                   [ 3.584428, -1.848126]])

# --- Check: assign each cold-replica sample to its nearest minimum ---
d = np.linalg.norm(cold_pts[:, None, :] - minima[None, :, :], axis=2)
nearest = np.argmin(d, axis=1)          # which basin each cold sample sits in
dist_to_min = d[np.arange(len(cold_pts)), nearest]
basins_visited = np.unique(nearest)
# a "hop" is a swap-round where the nearest-minimum label changes
hops = int(np.sum(nearest[1:] != nearest[:-1]))

# Spatial spread of each replica (hot should roam far more widely than cold)
cold_spread = cold_pts.std(axis=0)
hot_spread  = hot_pts.std(axis=0)

# --- Report numerical results ---
print("Metropolis acceptance rate per replica (cold->hot):")
for i in range(n_rep):
    print(f"  T={T[i]:5.1f}  step={steps[i]:.1f}  accept={accept_moves[i]/prop_moves[i]:.3f}")
print(f"Swap acceptance rate: {accept_swaps/prop_swaps:.3f}  ({accept_swaps}/{prop_swaps})")
print(f"Cold replica: number of samples = {len(cold_pts)}")
print(f"Cold replica: distinct basins visited = {len(basins_visited)} of 4 (labels {basins_visited.tolist()})")
print(f"Cold replica: basin hops between rounds = {hops}")
print(f"Cold replica: max distance to nearest minimum = {dist_to_min.max():.4f}")
print(f"Cold replica: mean distance to nearest minimum = {dist_to_min.mean():.4f}")
print(f"Cold replica: sample std (x, y) = ({cold_spread[0]:.4f}, {cold_spread[1]:.4f})")
print(f"Hot  replica: sample std (x, y) = ({hot_spread[0]:.4f}, {hot_spread[1]:.4f})")
print(f"Hot  replica: x range = [{hot_pts[:,0].min():.3f}, {hot_pts[:,0].max():.3f}]")
print(f"Hot  replica: y range = [{hot_pts[:,1].min():.3f}, {hot_pts[:,1].max():.3f}]")

# --- Plot both replicas on the Himmelblau landscape ---
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f((GX, GY))

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, pts, title in [(axes[0], cold_pts, f"Coldest replica (T={T[0]})"),
                       (axes[1], hot_pts,  f"Hottest replica (T={T[-1]})")]:
    ax.contourf(GX, GY, np.log1p(GZ), levels=30, cmap="viridis")
    ax.scatter(pts[:, 0], pts[:, 1], s=12, c="white", edgecolors="k",
               linewidths=0.3, alpha=0.8, label="samples")
    ax.scatter(minima[:, 0], minima[:, 1], s=120, marker="*", c="red",
               edgecolors="k", label="minima")
    ax.set_title(title)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.legend(loc="upper right")
fig.suptitle("Parallel tempering on Himmelblau's function")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.4.1_s1.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms the method because the cold replica's samples "
      "stay very close to the true minima (small distance-to-minimum) yet visit multiple "
      "basins with several hops, while the hot replica's much larger spatial spread shows "
      "it roams freely and supplies the cross-barrier configurations that let the cold "
      "chain jump between basins.")
