import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle switch in nondimensional form (Hill coefficient n=3, unit decay):
#   dx/dt = g/(1 + y^3) - x
#   dy/dt = h/(1 + x^3) - y
# X and Y mutually repress each other.
# ----------------------------------------------------------------------

g = 5.0
h = 5.0

def f(state):
    """Right-hand side of the ODE system. Returns [dx/dt, dy/dt]."""
    x, y = state
    dx = g / (1.0 + y**3) - x
    dy = h / (1.0 + x**3) - y
    return np.array([dx, dy])

def rk4_step(state, dt):
    """One explicit classical Runge-Kutta (RK4) step, written out by hand."""
    k1 = f(state)                    # slope at the start of the interval
    k2 = f(state + 0.5 * dt * k1)    # slope at the midpoint, using k1
    k3 = f(state + 0.5 * dt * k2)    # slope at the midpoint, using k2
    k4 = f(state + dt * k3)          # slope at the end, using k3
    # weighted average of the four slopes (1,2,2,1)/6
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

def integrate(state0, dt, n_steps):
    """Integrate the system forward in time and record the trajectory."""
    traj = np.empty((n_steps + 1, 2))
    traj[0] = state0
    s = np.array(state0, dtype=float)
    for i in range(n_steps):
        s = rk4_step(s, dt)
        traj[i + 1] = s
    return traj

# ----------------------------------------------------------------------
# Integrate from several initial conditions.
# ----------------------------------------------------------------------
dt = 0.01
n_steps = 3000  # t goes from 0 to 30, plenty for unit decay to settle

initial_conditions = [
    (4.0, 0.5),
    (3.0, 1.0),
    (2.0, 0.2),
    (0.5, 4.0),
    (1.0, 3.0),
    (0.2, 2.0),
    (1.0, 1.0),   # near the unstable symmetric point
    (2.0, 2.0),
    (0.1, 0.1),
    (4.5, 4.5),
]

trajectories = [integrate(ic, dt, n_steps) for ic in initial_conditions]

print("Final states reached from each initial condition:")
for ic, traj in zip(initial_conditions, trajectories):
    xf, yf = traj[-1]
    which = "x-high/y-low" if xf > yf else "x-low/y-high"
    print(f"  IC (x0={ic[0]:.2f}, y0={ic[1]:.2f}) -> (x={xf:.4f}, y={yf:.4f})  [{which}]")

# ----------------------------------------------------------------------
# CHECK: confirm two stable steady states by seeding trajectories in each
# basin and checking Jacobian eigenvalues at the converged fixed points.
# A fixed point is stable if all eigenvalues of the Jacobian have
# negative real part.
# ----------------------------------------------------------------------
def jacobian(state):
    """Analytic Jacobian of f at the given state."""
    x, y = state
    # d(dx/dt)/dx = -1 ; d(dx/dt)/dy = -3*g*y^2/(1+y^3)^2
    # d(dy/dt)/dx = -3*h*x^2/(1+x^3)^2 ; d(dy/dt)/dy = -1
    return np.array([
        [-1.0, -3.0 * g * y**2 / (1.0 + y**3)**2],
        [-3.0 * h * x**2 / (1.0 + x**3)**2, -1.0],
    ])

# Two long integrations, one seeded toward each basin, to land on fixed points.
fp_A = integrate((4.0, 0.1), dt, 6000)[-1]   # expect x-high/y-low
fp_B = integrate((0.1, 4.0), dt, 6000)[-1]   # expect x-low/y-high

print()
print("Steady state A (seeded x-high/y-low):")
print(f"  x* = {fp_A[0]:.6f}, y* = {fp_A[1]:.6f}")
print(f"  residual |f| = {np.linalg.norm(f(fp_A)):.3e}")
eigA = np.linalg.eigvals(jacobian(fp_A))
print(f"  Jacobian eigenvalues = {eigA[0]:.4f}, {eigA[1]:.4f}")
print(f"  stable? {np.all(eigA.real < 0)}")

print()
print("Steady state B (seeded x-low/y-high):")
print(f"  x* = {fp_B[0]:.6f}, y* = {fp_B[1]:.6f}")
print(f"  residual |f| = {np.linalg.norm(f(fp_B)):.3e}")
eigB = np.linalg.eigvals(jacobian(fp_B))
print(f"  Jacobian eigenvalues = {eigB[0]:.4f}, {eigB[1]:.4f}")
print(f"  stable? {np.all(eigB.real < 0)}")

print()
print("Two distinct stable steady states found:",
      not np.allclose(fp_A, fp_B) and np.all(eigA.real < 0) and np.all(eigB.real < 0))

# ----------------------------------------------------------------------
# Phase-plane plot of the trajectories.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))

for ic, traj in zip(initial_conditions, trajectories):
    xf = traj[-1, 0]
    yf = traj[-1, 1]
    color = "C0" if xf > yf else "C3"
    ax.plot(traj[:, 0], traj[:, 1], color=color, alpha=0.7, lw=1.2)
    ax.plot(traj[0, 0], traj[0, 1], "o", color=color, ms=5)

ax.plot(fp_A[0], fp_A[1], "k*", ms=18, label="stable SS x-high/y-low")
ax.plot(fp_B[0], fp_B[1], "k*", ms=18)
ax.plot(fp_B[0], fp_B[1], "P", color="gold", ms=10, label="stable SS x-low/y-high")
ax.plot(fp_A[0], fp_A[1], "P", color="cyan", ms=10)

ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title(f"Toggle switch phase plane (g={g}, h={h})\ntrajectories -> two stable steady states")
ax.legend(loc="upper right")
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.4.1_s4.png")

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ----------------------------------------------------------------------
print()
print("Explanation: Because two initial conditions converge to two different "
      "fixed points, each with all Jacobian eigenvalues having negative real "
      "part (hence locally stable), the system is confirmed bistable and the "
      "final state depends on which basin of attraction the initial condition lies in.")
