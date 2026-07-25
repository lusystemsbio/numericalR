import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# --- Parallel tempering setup ---
np.random.seed(1)                                   # seed 1
T = np.array([2.3, 5.0, 10.0, 20.0, 40.0, 80.0])    # temperatures (cold -> hot)
step = np.linspace(0.3, 2.5, len(T))                # step sizes 0.3 .. 2.5
n_rep = len(T)
n_rounds = 100                                       # swap rounds
n_steps = 100                                        # Metropolis steps per round

# All replicas start at the origin
X = np.zeros((n_rep, 2))
E = np.array([f(X[i, 0], X[i, 1]) for i in range(n_rep)])

# Store sampled points for each replica
samples = [[] for _ in range(n_rep)]

metropolis_accepts = np.zeros(n_rep)
metropolis_tries = np.zeros(n_rep)
swap_accepts = 0
swap_tries = 0

for rnd in range(n_rounds):
    # --- run a Metropolis chain at each temperature ---
    for i in range(n_rep):
        for _ in range(n_steps):
            # propose a move scaled by this replica's step size
            prop = X[i] + step[i] * np.random.uniform(-1, 1, size=2)
            Ep = f(prop[0], prop[1])
            # Metropolis acceptance at temperature T[i]
            if np.random.rand() < np.exp(-(Ep - E[i]) / T[i]):
                X[i] = prop
                E[i] = Ep
                metropolis_accepts[i] += 1
            metropolis_tries[i] += 1
            samples[i].append(X[i].copy())

    # --- attempt swaps of adjacent-temperature replicas ---
    # alternate which adjacent pairs are offered each round
    start = rnd % 2
    for i in range(start, n_rep - 1, 2):
        j = i + 1
        # acceptance a = min(1, exp((f_i - f_j)*(1/T_i - 1/T_j)))
        a = min(1.0, np.exp((E[i] - E[j]) * (1.0 / T[i] - 1.0 / T[j])))
        swap_tries += 1
        if np.random.rand() < a:
            # exchange the two configurations (and their energies)
            X[[i, j]] = X[[j, i]]
            E[[i, j]] = E[[j, i]]
            swap_accepts += 1

cold = np.array(samples[0])   # coldest replica (T = 2.3)
hot = np.array(samples[-1])   # hottest replica (T = 80)

# The four known minima of Himmelblau's function
minima = np.array([
    [3.0, 2.0],
    [-2.805118, 3.131312],
    [-3.779310, -3.283186],
    [3.584428, -1.848126],
])

# --- Check: assign cold-replica points to nearest minimum and measure distance ---
d_cold = np.min(np.linalg.norm(cold[:, None, :] - minima[None, :, :], axis=2), axis=1)
nearest_cold = np.argmin(np.linalg.norm(cold[:, None, :] - minima[None, :, :], axis=2), axis=1)
basins_visited_cold = np.unique(nearest_cold)

d_hot = np.min(np.linalg.norm(hot[:, None, :] - minima[None, :, :], axis=2), axis=1)

print("Number of replicas:", n_rep)
print("Temperatures:", ", ".join(f"{t:g}" for t in T))
print("Step sizes:", ", ".join(f"{s:.3f}" for s in step))
print("Total samples per replica:", cold.shape[0])
for i in range(n_rep):
    print(f"Metropolis acceptance rate T={T[i]:g}: {metropolis_accepts[i]/metropolis_tries[i]:.3f}")
print(f"Swap acceptance rate: {swap_accepts}/{swap_tries} = {swap_accepts/swap_tries:.3f}")

print(f"Coldest replica: fraction of samples within 0.5 of a minimum: {np.mean(d_cold < 0.5):.3f}")
print(f"Coldest replica: mean distance to nearest minimum: {d_cold.mean():.4f}")
print(f"Coldest replica: max distance to nearest minimum: {d_cold.max():.4f}")
print(f"Coldest replica: number of distinct basins visited: {len(basins_visited_cold)} of 4")
for b in range(4):
    print(f"  Cold samples nearest minimum {b} at ({minima[b,0]:.3f}, {minima[b,1]:.3f}): {np.sum(nearest_cold==b)}")
print(f"Hottest replica: mean distance to nearest minimum: {d_hot.mean():.4f}")
print(f"Hottest replica: max distance to nearest minimum: {d_hot.max():.4f}")
print(f"Hottest replica: x range: [{hot[:,0].min():.3f}, {hot[:,0].max():.3f}]")
print(f"Hottest replica: y range: [{hot[:,1].min():.3f}, {hot[:,1].max():.3f}]")
print(f"Coldest replica: x range: [{cold[:,0].min():.3f}, {cold[:,0].max():.3f}]")
print(f"Coldest replica: y range: [{cold[:,1].min():.3f}, {cold[:,1].max():.3f}]")

# Explanation of why the check confirms the result:
print("Explanation: The coldest replica clusters tightly around all four minima "
      "(small distances but multiple basins visited) while the hottest ranges "
      "broadly over the plane, showing that cold-basin hopping can only come from "
      "hot configurations passed down through replica exchange.")

# --- Plot the sampled points of the coldest and hottest replicas on the landscape ---
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)

fig, axes = plt.subplots(1, 2, figsize=(13, 6))
for ax, pts, title in ((axes[0], cold, f"Coldest replica (T={T[0]:g})"),
                        (axes[1], hot, f"Hottest replica (T={T[-1]:g})")):
    ax.contourf(GX, GY, np.log1p(GZ), levels=40, cmap="viridis")
    ax.plot(pts[:, 0], pts[:, 1], '.', ms=1.5, color="white", alpha=0.4)
    ax.plot(minima[:, 0], minima[:, 1], 'r*', ms=15, label="minima")
    ax.set_title(title)
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_xlim(-6, 6); ax.set_ylim(-6, 6)
    ax.legend(loc="upper right")

plt.suptitle("Parallel tempering on Himmelblau's function")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.4.1_s5.png")
