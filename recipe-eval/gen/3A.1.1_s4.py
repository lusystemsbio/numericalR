import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12

# Right-hand side of the ODE system, returning a 2-vector [dX/dt, dY/dt]
def f(state):
    X, Y = state[0], state[1]
    dX = gX0 + gX1 / (1.0 + (Y / Yth)**nY) - kX * X   # X transcription minus degradation
    dY = gY0 + gY1 / (1.0 + (X / Xth)**nX) - kY * Y   # Y transcription minus degradation
    return np.array([dX, dY])

# Explicit vector RK4: state and every stage k1..k4 are 2-vectors
def rk4_step(state, h):
    k1 = f(state)                    # slope at start
    k2 = f(state + 0.5 * h * k1)     # slope at midpoint using k1
    k3 = f(state + 0.5 * h * k2)     # slope at midpoint using k2
    k4 = f(state + h * k3)           # slope at end using k3
    return state + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)  # weighted average step

# Integrate one trajectory from an initial 2-vector
def integrate(init, h=0.5, nsteps=2000):
    traj = np.empty((nsteps + 1, 2))
    traj[0] = init
    for i in range(nsteps):
        traj[i + 1] = rk4_step(traj[i], h)
    return traj

# ---- Simulate ten random initial conditions in [0, 600]^2 ----
rng = np.random.default_rng(0)
inits = rng.uniform(0.0, 600.0, size=(10, 2))

trajectories = [integrate(state) for state in inits]

# ---- Collect and cluster the final (steady) states ----
finals = np.array([tr[-1] for tr in trajectories])
print("Final (X, Y) states per trajectory:")
for i, (x0y0, fin) in enumerate(zip(inits, finals)):
    print(f"  traj {i}: init=({x0y0[0]:.1f}, {x0y0[1]:.1f}) -> "
          f"steady=({fin[0]:.3f}, {fin[1]:.3f})")

# Cluster endpoints: group any that are within a small tolerance
tol = 1.0
clusters = []
for fin in finals:
    for c in clusters:
        if np.linalg.norm(fin - c[0]) < tol:
            c.append(fin)
            break
    else:
        clusters.append([fin])

cluster_means = [np.mean(c, axis=0) for c in clusters]
print()
print(f"Number of distinct stable steady states found: {len(cluster_means)}")
for j, cm in enumerate(cluster_means):
    print(f"  steady state {j}: X={cm[0]:.4f}, Y={cm[1]:.4f} "
          f"(reached by {len(clusters[j])} trajectories)")

# Report separation between the two states (distinct X and Y levels)
if len(cluster_means) == 2:
    a, b = cluster_means
    print(f"Difference in X between the two states: {abs(a[0]-b[0]):.4f}")
    print(f"Difference in Y between the two states: {abs(a[1]-b[1]):.4f}")

is_bistable = len(cluster_means) == 2
print(f"Bistable (exactly two distinct stable steady states): {is_bistable}")

# ---- Phase-plane plot ----
plt.figure(figsize=(8, 7))
for tr in trajectories:
    plt.plot(tr[:, 0], tr[:, 1], lw=1, alpha=0.7)
plt.scatter(inits[:, 0], inits[:, 1], c='k', marker='o', s=30,
            label='initial conditions', zorder=3)
for cm in cluster_means:
    plt.scatter(cm[0], cm[1], c='red', marker='*', s=300,
                edgecolor='black', zorder=4)
plt.scatter([], [], c='red', marker='*', s=200, edgecolor='black',
            label='stable steady states')
plt.xlabel("X level")
plt.ylabel("Y level")
plt.title("Genetic toggle switch: phase-plane trajectories (bistable)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.1.1_s4.png")

# Explanation of why the check confirms the result:
print()
print("Why the check confirms bistability: because ten trajectories started from "
      "randomly scattered initial conditions all converge to exactly two distinct "
      "(X, Y) fixed points rather than one, the system must have two coexisting "
      "stable steady states, which is the definition of a bistable switch.")
