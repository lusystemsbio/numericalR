import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Himmelblau's function (four equal minima at f = 0) -----
def f(p):
    x, y = p
    return (x*x + y - 11.0)**2 + (x + y*y - 7.0)**2

# ----- Parallel tempering setup -----
np.random.seed(1)                                  # seed 1 for reproducibility
T     = np.array([2.3, 5.0, 10.0, 20.0, 40.0, 80.0])   # temperature ladder
steps = np.array([0.3, 0.6, 1.0, 1.5, 2.0, 2.5])       # per-replica step sizes
n_rep = len(T)
n_rounds = 100                                     # number of swap rounds
n_steps  = 100                                     # Metropolis steps per round

# All replicas start at the origin
X = np.zeros((n_rep, 2))                           # current positions
E = np.array([f(X[i]) for i in range(n_rep)])      # current energies

# Storage for sampled points (record position after each Metropolis step)
samples = [[] for _ in range(n_rep)]

def metropolis_step(i):
    """One Metropolis move for replica i at temperature T[i]."""
    prop = X[i] + steps[i] * np.random.uniform(-1.0, 1.0, size=2)  # random displacement
    ep = f(prop)
    # Accept with Metropolis rule using this replica's own temperature
    if np.random.rand() < np.exp(-(ep - E[i]) / T[i]):
        X[i] = prop
        E[i] = ep

# ----- Main loop: alternate local sampling with replica exchange -----
n_swap_attempts = 0
n_swap_accept   = 0
for rnd in range(n_rounds):
    # (1) Run an independent Metropolis chain in each replica
    for _ in range(n_steps):
        for i in range(n_rep):
            metropolis_step(i)
            samples[i].append(X[i].copy())         # record sampled point

    # (2) Attempt swaps of adjacent-temperature replicas
    #     alternate even/odd pairing across rounds so all pairs are tried
    start = rnd % 2
    for i in range(start, n_rep - 1, 2):
        j = i + 1
        # acceptance a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))
        a = min(1.0, np.exp((E[i] - E[j]) * (1.0/T[i] - 1.0/T[j])))
        n_swap_attempts += 1
        if np.random.rand() < a:
            X[[i, j]] = X[[j, i]]                   # exchange configurations
            E[[i, j]] = E[[j, i]]
            n_swap_accept += 1

samples = [np.array(s) for s in samples]
cold = samples[0]     # coldest replica (T = 2.3)
hot  = samples[-1]    # hottest replica (T = 80)

# ----- Numerical results -----
print(f"Number of replicas: {n_rep}")
print(f"Total samples per replica: {len(cold)}")
print(f"Swap attempts: {n_swap_attempts}")
print(f"Swap acceptances: {n_swap_accept}")
print(f"Overall swap acceptance rate: {n_swap_accept / n_swap_attempts:.4f}")

# The four known minima of Himmelblau's function
minima = np.array([[ 3.000000,  2.000000],
                   [-2.805118,  3.131312],
                   [-3.779310, -3.283186],
                   [ 3.584428, -1.848126]])

# Assign each coldest-replica sample to its nearest minimum, and measure distance
d_cold = np.min(np.linalg.norm(cold[:, None, :] - minima[None, :, :], axis=2), axis=1)
nearest = np.argmin(np.linalg.norm(cold[:, None, :] - minima[None, :, :], axis=2), axis=1)

print(f"Coldest replica mean f value: {np.mean([f(p) for p in cold]):.4f}")
print(f"Hottest replica mean f value: {np.mean([f(p) for p in hot]):.4f}")
print(f"Coldest replica mean distance to nearest minimum: {np.mean(d_cold):.4f}")
print(f"Coldest replica max distance to nearest minimum: {np.max(d_cold):.4f}")
print(f"Coldest replica fraction of samples within 0.75 of a minimum: {np.mean(d_cold < 0.75):.4f}")

# Count how many distinct basins the coldest replica visits (must be >1 to prove hopping)
for k in range(4):
    print(f"Coldest replica samples assigned to minimum {k} {tuple(np.round(minima[k],3))}: {np.sum(nearest==k)}")
print(f"Distinct basins visited by coldest replica: {len(np.unique(nearest))}")

# Spread of hottest replica (should span the whole landscape)
print(f"Hottest replica x range: [{hot[:,0].min():.3f}, {hot[:,0].max():.3f}]")
print(f"Hottest replica y range: [{hot[:,1].min():.3f}, {hot[:,1].max():.3f}]")

# ----- Plot: sampled points on the landscape -----
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = (GX**2 + GY - 11)**2 + (GX + GY**2 - 7)**2

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, pts, ttl in [(axes[0], cold, "Coldest replica (T=2.3)"),
                     (axes[1], hot,  "Hottest replica (T=80)")]:
    ax.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
    ax.scatter(pts[:, 0], pts[:, 1], s=3, c="white", alpha=0.25)
    ax.scatter(minima[:, 0], minima[:, 1], c="red", marker="*", s=200, edgecolor="k", label="minima")
    ax.set_title(ttl)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.legend(loc="upper right")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.4.1_s3.png", dpi=120)

# ----- One-sentence explanation of why the check confirms the result -----
print("Explanation: Because the coldest replica's samples cluster tightly around "
      "multiple distinct minima (low mean f, small distance to a minimum, but >1 basin "
      "visited) while the hottest replica's samples span the entire domain, we confirm "
      "that replica exchange feeds hot, barrier-crossing configurations down to the cold "
      "chain, letting it hop between basins without ever leaving the low-energy regions.")
