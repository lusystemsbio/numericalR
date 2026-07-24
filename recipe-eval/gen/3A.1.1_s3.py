import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Genetic toggle switch: genes X and Y mutually repress each other.
# State is a 2-vector u = [X, Y].  We integrate with an explicit
# vector-form RK4 (each stage k1..k4 is itself a 2-vector).
# ----------------------------------------------------------------------

# Model parameters
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12


def f(u):
    """Right-hand side of the ODE system; returns a 2-vector du/dt."""
    X, Y = u[0], u[1]
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X   # X transcription minus decay
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y   # Y transcription minus decay
    return np.array([dX, dY])


def rk4_step(u, dt):
    """One explicit RK4 step, carrying the state and every stage as 2-vectors."""
    k1 = f(u)                    # slope at start
    k2 = f(u + 0.5 * dt * k1)    # slope at midpoint using k1
    k3 = f(u + 0.5 * dt * k2)    # slope at midpoint using k2
    k4 = f(u + dt * k3)          # slope at end using k3
    return u + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)  # weighted average


def integrate(u0, dt, n_steps):
    """Integrate from initial 2-vector u0, returning the full trajectory."""
    traj = np.empty((n_steps + 1, 2))
    traj[0] = u0
    for i in range(n_steps):
        traj[i + 1] = rk4_step(traj[i], dt)
    return traj


# ----------------------------------------------------------------------
# Simulate ten random initial conditions in [0, 600]^2
# ----------------------------------------------------------------------
rng = np.random.default_rng(0)
n_ic = 10
dt = 0.5
n_steps = 2000  # long enough to reach steady state (T = 1000)

inits = rng.uniform(0.0, 600.0, size=(n_ic, 2))
trajectories = []
final_states = []

for j in range(n_ic):
    tr = integrate(inits[j], dt, n_steps)
    trajectories.append(tr)
    final_states.append(tr[-1])
    print(f"IC {j:2d}: start (X0={inits[j,0]:7.2f}, Y0={inits[j,1]:7.2f})  ->  "
          f"final (X={tr[-1,0]:8.3f}, Y={tr[-1,1]:8.3f})")

final_states = np.array(final_states)

# ----------------------------------------------------------------------
# Check: cluster the final states and confirm exactly two distinct stable
# steady states (different X and Y levels).
# ----------------------------------------------------------------------
tol = 1e-2  # states within this distance are the "same" steady state
clusters = []          # list of representative steady-state 2-vectors
labels = np.empty(n_ic, dtype=int)

for j in range(n_ic):
    assigned = False
    for c, rep in enumerate(clusters):
        if np.linalg.norm(final_states[j] - rep) < tol:
            labels[j] = c
            assigned = True
            break
    if not assigned:
        clusters.append(final_states[j])
        labels[j] = len(clusters) - 1

clusters = np.array(clusters)
n_stable = len(clusters)

print()
print(f"Number of distinct stable steady states found: {n_stable}")
for c, rep in enumerate(clusters):
    count = int(np.sum(labels == c))
    print(f"Steady state {c}: X={rep[0]:8.3f}, Y={rep[1]:8.3f}  "
          f"(reached by {count} of {n_ic} trajectories)")

if n_stable == 2:
    sepX = abs(clusters[0, 0] - clusters[1, 0])
    sepY = abs(clusters[0, 1] - clusters[1, 1])
    print(f"Separation between the two states: dX={sepX:.3f}, dY={sepY:.3f}")
    bistable = sepX > tol and sepY > tol
    print(f"Two states differ in BOTH X and Y (distinct levels): {bistable}")
    print(f"System is bistable: {bistable}")
else:
    print("System is bistable: False")

# One-sentence explanation of why this check confirms bistability:
print()
print("Explanation: Because every one of the many random initial conditions "
      "settles onto exactly one of two distinct stable fixed points (differing "
      "in both X and Y), the system has two coexisting stable states and is "
      "therefore bistable.")

# ----------------------------------------------------------------------
# Phase-plane plot of the ten trajectories converging to the two states
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 7))
colors = ['C0', 'C3']  # color trajectories by which steady state they reach
for j in range(n_ic):
    c = labels[j] if labels[j] < 2 else 0
    tr = trajectories[j]
    plt.plot(tr[:, 0], tr[:, 1], color=colors[c], alpha=0.6, lw=1.0)
    plt.plot(tr[0, 0], tr[0, 1], 'o', color=colors[c], ms=5)  # start

# Mark the stable steady states
for c, rep in enumerate(clusters):
    plt.plot(rep[0], rep[1], 'k*', ms=18,
             label=f"Stable state {c} (X={rep[0]:.1f}, Y={rep[1]:.1f})")

plt.xlabel("X")
plt.ylabel("Y")
plt.title("Genetic toggle switch: 10 trajectories converging to 2 stable states")
plt.legend(loc="best")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.1.1_s3.png")
print()
print("Saved phase-plane plot to 3A.1.1_s3.png")
