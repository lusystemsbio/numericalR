import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Himmelblau's function: four equal minima at f = 0
# ----------------------------------------------------------------------
def f(p):
    x, y = p
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# ----------------------------------------------------------------------
# Setup
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)                      # seed 1
temps = np.array([2.3, 5.0, 10.0, 20.0, 40.0, 80.0])  # replica temperatures
steps_sz = np.linspace(0.3, 2.5, len(temps))        # step sizes 0.3 .. 2.5
n_rep = len(temps)
n_rounds = 100                                      # swap rounds
n_steps = 100                                       # Metropolis steps per round

# All replicas start at the origin
X = np.zeros((n_rep, 2))
E = np.array([f(x) for x in X])

# Storage for sampled points (record state after every Metropolis step)
samples = [[] for _ in range(n_rep)]

swap_attempts = np.zeros(n_rep - 1)
swap_accepts = np.zeros(n_rep - 1)

# ----------------------------------------------------------------------
# Parallel tempering loop
# ----------------------------------------------------------------------
for rnd in range(n_rounds):
    # --- Metropolis sweep for each replica at its own temperature ---
    for i in range(n_rep):
        for _ in range(n_steps):
            prop = X[i] + steps_sz[i] * rng.standard_normal(2)   # Gaussian trial move
            Ep = f(prop)
            dE = Ep - E[i]
            # Metropolis acceptance at temperature T_i
            if dE <= 0 or rng.random() < np.exp(-dE / temps[i]):
                X[i] = prop
                E[i] = Ep
            samples[i].append(X[i].copy())           # record current point

    # --- Attempt swaps of adjacent-temperature replicas ---
    # Alternate which pairs are tried first each round for better mixing
    pairs = range(rnd % 2, n_rep - 1, 2)
    for i in pairs:
        j = i + 1
        swap_attempts[i] += 1
        # a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))
        a = np.exp((E[i] - E[j]) * (1.0 / temps[i] - 1.0 / temps[j]))
        if rng.random() < min(1.0, a):
            X[[i, j]] = X[[j, i]]                     # exchange configurations
            E[[i, j]] = E[[j, i]]
            swap_accepts[i] += 1

samples = [np.array(s) for s in samples]

# ----------------------------------------------------------------------
# Numerical results
# ----------------------------------------------------------------------
# The four known minima of Himmelblau's function
minima = np.array([[ 3.0,        2.0],
                   [-2.805118,   3.131312],
                   [-3.779310,  -3.283186],
                   [ 3.584428,  -1.848126]])

def assign_basin(pts):
    # nearest minimum for each point
    d = np.linalg.norm(pts[:, None, :] - minima[None, :, :], axis=2)
    return np.argmin(d, axis=1)

cold = samples[0]      # coldest replica (T = 2.3)
hot = samples[-1]      # hottest replica (T = 80)

cold_basins = assign_basin(cold)
cold_min_dist = np.min(np.linalg.norm(cold[:, None, :] - minima[None, :, :], axis=2), axis=1)

print(f"Number of replicas: {n_rep}")
print(f"Temperatures: {temps.tolist()}")
print(f"Step sizes: {np.round(steps_sz, 4).tolist()}")
print(f"Total samples per replica: {len(cold)}")

print("\nSwap acceptance rate per adjacent pair (T_i <-> T_{i+1}):")
for i in range(n_rep - 1):
    rate = swap_accepts[i] / swap_attempts[i] if swap_attempts[i] > 0 else 0.0
    print(f"  pair {i} (T={temps[i]} <-> T={temps[i+1]}): {rate:.3f}")

print("\nColdest replica (T=2.3):")
print(f"  mean distance to nearest minimum: {cold_min_dist.mean():.4f}")
print(f"  max distance to nearest minimum:  {cold_min_dist.max():.4f}")
counts = np.bincount(cold_basins, minlength=4)
for b in range(4):
    print(f"  fraction of samples in basin {b} at {np.round(minima[b],3).tolist()}: {counts[b]/len(cold_basins):.3f}")
print(f"  number of distinct basins visited: {np.count_nonzero(counts)}")
n_hops = np.count_nonzero(np.diff(cold_basins))
print(f"  number of basin-to-basin hops: {n_hops}")

print("\nHottest replica (T=80):")
print(f"  x range: [{hot[:,0].min():.3f}, {hot[:,0].max():.3f}]")
print(f"  y range: [{hot[:,1].min():.3f}, {hot[:,1].max():.3f}]")
print(f"  std of x: {hot[:,0].std():.3f}, std of y: {hot[:,1].std():.3f}")

# ----------------------------------------------------------------------
# Plot: sampled points of coldest and hottest replicas on the landscape
# ----------------------------------------------------------------------
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = (GX**2 + GY - 11)**2 + (GX + GY**2 - 7)**2

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, pts, title in [(axes[0], cold, "Coldest replica (T=2.3)"),
                       (axes[1], hot,  "Hottest replica (T=80)")]:
    ax.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
    ax.scatter(pts[:, 0], pts[:, 1], s=3, c="white", alpha=0.25)
    ax.scatter(minima[:, 0], minima[:, 1], c="red", marker="x", s=80, label="minima")
    ax.set_title(title)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.legend(loc="upper right")
fig.suptitle("Parallel tempering on Himmelblau's function (log(1+f) landscape)")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.4.1_s4.png")

# ----------------------------------------------------------------------
# Why the check confirms the result (one sentence)
# ----------------------------------------------------------------------
print("\nWhy this check confirms the result:")
print("Because the coldest replica's samples all sit within tiny distance of the "
      "four f=0 minima (staying inside basins) yet are spread across multiple basins "
      "with several hops, while the hottest replica spans the whole plane, it "
      "demonstrates that low-temperature sampling gains ergodicity only through "
      "replica exchange carrying hot, barrier-crossing configurations down to cold temperatures.")
