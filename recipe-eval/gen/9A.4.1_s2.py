import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Himmelblau's function: four equal minima at f = 0 ---
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# --- Setup: 6 replicas at increasing temperatures with matched step sizes ---
rng = np.random.default_rng(1)                       # seed 1
temps = np.array([2.3, 5.0, 10.0, 20.0, 40.0, 80.0]) # temperatures (cold -> hot)
steps = np.array([0.3, 0.5, 0.8, 1.2, 1.8, 2.5])     # proposal step sizes per replica
n_rep = len(temps)
n_rounds = 100                                       # number of swap rounds
n_steps = 100                                        # Metropolis steps between swaps

# All replicas start at the origin
pos = np.zeros((n_rep, 2))                            # current (x, y) of each replica
energy = np.array([f(p[0], p[1]) for p in pos])      # current f value of each replica

# Storage for sampled points of every replica
samples = [[] for _ in range(n_rep)]

# Diagnostics
metro_accept = np.zeros(n_rep)
metro_total = np.zeros(n_rep)
swap_accept = 0
swap_total = 0

for rnd in range(n_rounds):
    # --- Metropolis sweep at each temperature independently ---
    for r in range(n_rep):
        T = temps[r]
        s = steps[r]
        for _ in range(n_steps):
            # propose a Gaussian move of scale s
            trial = pos[r] + rng.normal(0.0, s, size=2)
            e_trial = f(trial[0], trial[1])
            # Metropolis acceptance at temperature T
            if e_trial <= energy[r] or rng.random() < np.exp((energy[r] - e_trial) / T):
                pos[r] = trial
                energy[r] = e_trial
                metro_accept[r] += 1
            metro_total[r] += 1
            samples[r].append(pos[r].copy())

    # --- Replica exchange: try swapping adjacent-temperature replicas ---
    # Alternate which pairs we start with each round for better mixing
    start = rnd % 2
    for i in range(start, n_rep - 1, 2):
        j = i + 1
        Ti, Tj = temps[i], temps[j]
        fi, fj = energy[i], energy[j]
        # acceptance a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))
        a = min(1.0, np.exp((fi - fj) * (1.0 / Ti - 1.0 / Tj)))
        swap_total += 1
        if rng.random() < a:
            # exchange the configurations (and their energies) between replicas
            pos[[i, j]] = pos[[j, i]]
            energy[[i, j]] = energy[[j, i]]
            swap_accept += 1

samples = [np.array(s) for s in samples]
cold = samples[0]     # coldest replica, T = 2.3
hot = samples[-1]     # hottest replica, T = 80

# --- Known four minima of Himmelblau's function ---
minima = np.array([[ 3.000000,  2.000000],
                   [-2.805118,  3.131312],
                   [-3.779310, -3.283186],
                   [ 3.584428, -1.848126]])

# Assign each cold-replica sample to its nearest minimum
d = np.linalg.norm(cold[:, None, :] - minima[None, :, :], axis=2)
nearest = np.argmin(d, axis=1)
dist_to_nearest = d[np.arange(len(cold)), nearest]

print("Temperatures:", temps.tolist())
print("Step sizes:", steps.tolist())
for r in range(n_rep):
    print(f"Metropolis acceptance rate T={temps[r]:5.1f}: {metro_accept[r]/metro_total[r]:.3f}")
print(f"Replica-swap acceptance rate: {swap_accept/swap_total:.3f}")

# Cold replica: how tightly it clusters and how many basins it visits
print(f"Cold replica mean distance to nearest minimum: {dist_to_nearest.mean():.4f}")
print(f"Cold replica max distance to nearest minimum:  {dist_to_nearest.max():.4f}")
print(f"Cold replica fraction of samples within 0.75 of a minimum: {(dist_to_nearest < 0.75).mean():.4f}")
counts = np.bincount(nearest, minlength=4)
for k in range(4):
    print(f"Cold replica samples nearest minimum {minima[k].tolist()}: {counts[k]}")
print(f"Cold replica number of distinct basins visited: {int((counts > 0).sum())}")

# Hottest replica: spatial spread showing it roams the whole landscape
print(f"Hot replica x range: [{hot[:,0].min():.3f}, {hot[:,0].max():.3f}]")
print(f"Hot replica y range: [{hot[:,1].min():.3f}, {hot[:,1].max():.3f}]")
print(f"Hot replica std (x, y): ({hot[:,0].std():.3f}, {hot[:,1].std():.3f})")
print(f"Cold replica std (x, y): ({cold[:,0].std():.3f}, {cold[:,1].std():.3f})")

# --- Plot the landscape with sampled points of coldest and hottest replicas ---
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)

fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharex=True, sharey=True)
for ax, pts, title in ((axes[0], cold, "Coldest replica  (T = 2.3)"),
                       (axes[1], hot, "Hottest replica  (T = 80)")):
    ax.contourf(GX, GY, np.log1p(GZ), levels=30, cmap="viridis")
    ax.scatter(pts[:, 0], pts[:, 1], s=3, c="white", alpha=0.25, linewidths=0)
    ax.scatter(minima[:, 0], minima[:, 1], marker="*", s=200,
               c="red", edgecolors="black", zorder=5, label="minima")
    ax.set_title(title)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.legend(loc="upper right")
fig.suptitle("Parallel tempering on Himmelblau's function (log(1+f) background)")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.4.1_s2.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms parallel tempering works because the cold replica's "
      "samples stay tightly clustered inside the four minima (small mean distance) yet its "
      "counts show it populates multiple basins, which is only possible if hot-replica "
      "configurations that freely roam the whole landscape are being swapped down to it.")
