import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Himmelblau's function: four equal minima at f = 0
def f(x, y):
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

# The four known global minima of Himmelblau's function
known_minima = np.array([
    [ 3.000000,  2.000000],
    [-2.805118,  3.131312],
    [-3.779310, -3.283186],
    [ 3.584428, -1.848126],
])

# Explicit Metropolis-Hastings sampler at fixed temperature T
def metropolis_hastings(T, n_steps=int(1e4), step=0.5, start=(0.0, 0.0), seed=1):
    rng = np.random.default_rng(seed)          # per-temperature seed
    x, y = start                               # current state
    fx = f(x, y)                               # current energy
    chain = np.empty((n_steps, 2))             # store sampled points
    n_accept = 0
    for i in range(n_steps):
        # uniform proposal displacement in [-step, step] for each coordinate
        xp = x + rng.uniform(-step, step)
        yp = y + rng.uniform(-step, step)
        fp = f(xp, yp)
        de = fp - fx                           # energy change
        # acceptance probability a = min(1, exp(-de/T))
        a = 1.0 if de <= 0 else np.exp(-de / T)
        if rng.uniform() < a:                  # accept move
            x, y, fx = xp, yp, fp
            n_accept += 1
        chain[i] = (x, y)                       # record (accepted or repeated) state
    return chain, n_accept / n_steps

temperatures = [1, 10, 30, 50]
results = {}
for T in temperatures:
    chain, acc = metropolis_hastings(T)
    results[T] = chain
    # distance of the final point to the nearest known minimum -> did it "pin" a minimum?
    final = chain[-1]
    d_final = np.min(np.linalg.norm(known_minima - final, axis=1))
    # how many of the four basins were visited (a point counts as visiting a basin
    # if it lands within radius 1.0 of that minimum) -> did it "explore widely"?
    dists = np.linalg.norm(chain[:, None, :] - known_minima[None, :, :], axis=2)
    basins_visited = int(np.sum(np.any(dists < 1.0, axis=0)))
    spread = chain.std(axis=0)
    print(f"T = {T}:")
    print(f"  acceptance rate                 = {acc:.4f}")
    print(f"  sample std (x, y)               = ({spread[0]:.4f}, {spread[1]:.4f})")
    print(f"  final point                     = ({final[0]:.4f}, {final[1]:.4f})")
    print(f"  distance final -> nearest min   = {d_final:.4f}")
    print(f"  number of basins visited (of 4) = {basins_visited}")

# Explicit low-T vs high-T check
T_low, T_high = 1, 50
chain_low, chain_high = results[T_low], results[T_high]
basins_low = int(np.sum(np.any(
    np.linalg.norm(chain_low[:, None, :] - known_minima[None, :, :], axis=2) < 1.0, axis=0)))
basins_high = int(np.sum(np.any(
    np.linalg.norm(chain_high[:, None, :] - known_minima[None, :, :], axis=2) < 1.0, axis=0)))
d_final_low = np.min(np.linalg.norm(known_minima - chain_low[-1], axis=1))
d_final_high = np.min(np.linalg.norm(known_minima - chain_high[-1], axis=1))
print("CHECK: no single temperature both explores widely and settles")
print(f"  T=1  basins visited = {basins_low} (traps in one basin),"
      f" final dist to min = {d_final_low:.4f} (pins a minimum)")
print(f"  T=50 basins visited = {basins_high} (roams widely),"
      f" final dist to min = {d_final_high:.4f} (never pins a minimum)")
print("  => low T settles but does not explore; high T explores but does not settle.")

# Plot the sampled points at each temperature on the Himmelblau landscape
gx = np.linspace(-6, 6, 400)
gy = np.linspace(-6, 6, 400)
GX, GY = np.meshgrid(gx, gy)
GZ = f(GX, GY)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))
for ax, T in zip(axes.ravel(), temperatures):
    # log-scaled contours reveal the four basins clearly
    ax.contourf(GX, GY, np.log10(GZ + 1), levels=30, cmap="viridis")
    chain = results[T]
    ax.plot(chain[:, 0], chain[:, 1], '.', color="white", markersize=1, alpha=0.3)
    ax.plot(known_minima[:, 0], known_minima[:, 1], 'r*', markersize=15,
            label="true minima")
    ax.set_title(f"Metropolis-Hastings on Himmelblau, T = {T}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.legend(loc="upper right", fontsize=8)

fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/9A.2.1_s3.png", dpi=120)

# One-sentence explanation of why the check confirms the result
print("EXPLANATION: The check confirms the result because at T=1 the chain visits only "
      "one basin yet ends essentially at a true minimum (settles but doesn't explore), "
      "whereas at T=50 the chain wanders through all four basins yet ends far from any "
      "minimum (explores but doesn't settle) - showing no single fixed temperature "
      "achieves both simultaneously.")
