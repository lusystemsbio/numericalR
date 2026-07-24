import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12

# ---- Vector field: takes 2-vector state, returns 2-vector derivative ----
def f(s):
    X, Y = s[0], s[1]
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# ---- Explicit vector RK4: state and each stage are 2-vectors ----
def rk4_step(s, h):
    k1 = f(s)              # slope at start
    k2 = f(s + 0.5 * h * k1)  # slope at midpoint using k1
    k3 = f(s + 0.5 * h * k2)  # slope at midpoint using k2
    k4 = f(s + h * k3)        # slope at end using k3
    return s + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # weighted average

def integrate(s0, h, nsteps):
    traj = np.empty((nsteps + 1, 2))
    traj[0] = s0
    s = s0.copy()
    for i in range(nsteps):   # march forward one RK4 step at a time
        s = rk4_step(s, h)
        traj[i + 1] = s
    return traj

# ---- Simulate ten random initial conditions in [0, 600] ----
rng = np.random.default_rng(0)
inits = rng.uniform(0.0, 600.0, size=(10, 2))
h, T = 0.5, 800.0
nsteps = int(T / h)

trajectories = []
endpoints = []
for j in range(10):
    tr = integrate(inits[j], h, nsteps)
    trajectories.append(tr)
    endpoints.append(tr[-1])
    print(f"IC {j}: start (X0={inits[j,0]:.2f}, Y0={inits[j,1]:.2f}) -> "
          f"steady (X={tr[-1,0]:.4f}, Y={tr[-1,1]:.4f})")

endpoints = np.array(endpoints)

# ---- Bistability check: cluster the endpoints into distinct steady states ----
# Group endpoints that are within a small tolerance of each other.
tol = 1.0
clusters = []
for p in endpoints:
    placed = False
    for c in clusters:
        if np.linalg.norm(p - c[0]) < tol:
            c.append(p)
            placed = True
            break
    if not placed:
        clusters.append([p])

centers = [np.mean(c, axis=0) for c in clusters]
print(f"\nNumber of distinct stable steady states found: {len(centers)}")
for i, (c, members) in enumerate(zip(centers, clusters)):
    print(f"Steady state {i+1}: X={c[0]:.4f}, Y={c[1]:.4f}  "
          f"(reached by {len(members)} of 10 trajectories)")

if len(centers) == 2:
    sep = np.linalg.norm(centers[0] - centers[1])
    print(f"Separation between the two steady states: {sep:.4f}")
    print(f"dX between states: {abs(centers[0][0]-centers[1][0]):.4f}")
    print(f"dY between states: {abs(centers[0][1]-centers[1][1]):.4f}")
    print("Bistable: True (every trajectory settled on one of exactly two distinct states)")
else:
    print("Bistable: False (did not find exactly two distinct states)")

# One-sentence explanation:
print("\nWhy this confirms the result: because all ten trajectories, started from "
      "widely scattered initial conditions, collapse onto exactly two well-separated "
      "fixed points (distinct X and Y levels), the system has two coexisting stable "
      "steady states, which is the definition of a bistable toggle switch.")

# ---- Phase-plane plot ----
plt.figure(figsize=(8, 7))
for j, tr in enumerate(trajectories):
    plt.plot(tr[:, 0], tr[:, 1], lw=1, alpha=0.8)
    plt.plot(tr[0, 0], tr[0, 1], 'o', color='gray', ms=5)  # start markers
for c in centers:
    plt.plot(c[0], c[1], 'k*', ms=20, zorder=5)  # stable steady states
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Genetic toggle switch: 10 trajectories converging to two stable states")
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.1.1_s2.png", dpi=120, bbox_inches="tight")
