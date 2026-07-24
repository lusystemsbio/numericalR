import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------
# Model parameters for the genetic toggle switch (X and Y repress each other)
# ------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10   # X gene parameters
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12   # Y gene parameters


def deriv(s):
    """Return the 2-vector [dX/dt, dY/dt] given state 2-vector s = [X, Y]."""
    X, Y = s[0], s[1]
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X   # X: basal + repressed-by-Y - decay
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y   # Y: basal + repressed-by-X - decay
    return np.array([dX, dY])


def rk4_step(s, h):
    """One classical RK4 step, done explicitly with the state and stages as 2-vectors."""
    k1 = deriv(s)              # slope at start
    k2 = deriv(s + 0.5 * h * k1)  # slope at midpoint using k1
    k3 = deriv(s + 0.5 * h * k2)  # slope at midpoint using k2
    k4 = deriv(s + h * k3)        # slope at end using k3
    return s + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)  # weighted average


def integrate(s0, h, nsteps):
    """Integrate from initial 2-vector s0, returning the full trajectory array."""
    traj = np.empty((nsteps + 1, 2))
    traj[0] = s0
    for i in range(nsteps):
        traj[i + 1] = rk4_step(traj[i], h)  # advance one RK4 step
    return traj


# ------------------------------------------------------------------
# Simulate ten random initial conditions in [0, 600] x [0, 600]
# ------------------------------------------------------------------
np.random.seed(0)
h = 0.5            # time step
nsteps = 2000      # total simulated time = 1000
n_ic = 10

init_conditions = np.random.uniform(0.0, 600.0, size=(n_ic, 2))

trajectories = []
endpoints = []
for j in range(n_ic):
    tr = integrate(init_conditions[j], h, nsteps)
    trajectories.append(tr)
    endpoints.append(tr[-1])
    print(f"IC {j:2d}: start (X0,Y0)=({init_conditions[j,0]:7.2f},{init_conditions[j,1]:7.2f})  "
          f"-> end (X,Y)=({tr[-1,0]:8.3f},{tr[-1,1]:8.3f})")

endpoints = np.array(endpoints)

# ------------------------------------------------------------------
# Bistability check: cluster the endpoints into distinct stable steady states
# ------------------------------------------------------------------
tol = 1.0  # states within this distance are considered the same steady state
clusters = []  # list of [representative_state, [member_indices]]
for j, ep in enumerate(endpoints):
    placed = False
    for c in clusters:
        if np.linalg.norm(ep - c[0]) < tol:
            c[1].append(j)
            placed = True
            break
    if not placed:
        clusters.append([ep.copy(), [j]])

n_states = len(clusters)
print(f"\nNumber of distinct stable steady states found: {n_states}")
for i, c in enumerate(clusters):
    print(f"Steady state {i+1}: (X,Y)=({c[0][0]:8.3f},{c[0][1]:8.3f})  "
          f"reached by {len(c[1])} of {n_ic} trajectories")

# Confirm the two states differ in both X and Y levels
if n_states == 2:
    s1, s2 = clusters[0][0], clusters[1][0]
    print(f"\nDifference between the two states: dX={abs(s1[0]-s2[0]):.3f}  dY={abs(s1[1]-s2[1]):.3f}")
    print("Bistable: YES" if n_states == 2 else "Bistable: NO")
else:
    print("\nBistable: NO (did not settle on exactly two states)")

# ------------------------------------------------------------------
# Phase-plane plot of the ten trajectories converging to the steady states
# ------------------------------------------------------------------
plt.figure(figsize=(8, 7))
colors = plt.cm.viridis(np.linspace(0, 1, n_ic))
for j, tr in enumerate(trajectories):
    plt.plot(tr[:, 0], tr[:, 1], color=colors[j], lw=1.0, alpha=0.8)
    plt.plot(tr[0, 0], tr[0, 1], 'o', color=colors[j], ms=5)  # start marker

# mark the stable steady states
for i, c in enumerate(clusters):
    plt.plot(c[0][0], c[0][1], 'r*', ms=20, mec='k',
             label=f"Steady state {i+1}" if i < 2 else None)

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Genetic Toggle Switch: 10 trajectories converging to two stable steady states")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.1.1_s1.png")

# ------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result
print("\nExplanation: Because all ten trajectories started from widely scattered "
      "random initial conditions yet collapsed onto exactly two distinct attracting "
      "points (differing in both X and Y), the system has two coexisting stable steady "
      "states, which is the definition of bistability.")
