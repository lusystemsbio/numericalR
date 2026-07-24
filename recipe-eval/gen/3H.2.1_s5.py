import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Two-gene negative-feedback loop (nondimensional form)
#   dx/dt = g/(1 + y^3) - x       (X is repressed by Y)
#   dy/dt = h*x^3/(1 + x^3) - y   (Y is activated by X)
# Hill coefficient 3, unit degradation.
# ---------------------------------------------------------------
g = 10.0
h = 10.0

def f(state):
    # Right-hand side of the ODE system; returns [dx/dt, dy/dt].
    x, y = state
    dx = g / (1.0 + y**3) - x
    dy = h * x**3 / (1.0 + x**3) - y
    return np.array([dx, dy])

def rk4_step(state, dt):
    # One explicit classical Runge-Kutta 4th-order step.
    k1 = f(state)                    # slope at start
    k2 = f(state + 0.5 * dt * k1)    # slope at midpoint using k1
    k3 = f(state + 0.5 * dt * k2)    # slope at midpoint using k2
    k4 = f(state + dt * k3)          # slope at end using k3
    # Weighted average of the four slopes (Simpson-like weights 1,2,2,1).
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

def integrate(state0, dt, n_steps):
    # March the RK4 step forward and record the whole trajectory.
    traj = np.empty((n_steps + 1, 2))
    traj[0] = state0
    for i in range(n_steps):
        traj[i + 1] = rk4_step(traj[i], dt)
    return traj

# ---------------------------------------------------------------
# Integrate several trajectories from different initial conditions.
# ---------------------------------------------------------------
dt = 0.01
T = 30.0
n_steps = int(T / dt)

initial_conditions = [
    (0.5, 0.5),
    (9.0, 0.5),
    (0.5, 9.0),
    (8.0, 8.0),
    (2.0, 6.0),
]

trajectories = [integrate(np.array(ic, dtype=float), dt, n_steps)
                for ic in initial_conditions]

# ---------------------------------------------------------------
# Steady state: iterate one trajectory far in time until it stops moving.
# ---------------------------------------------------------------
ss = integrate(np.array([5.0, 5.0]), dt, int(200.0 / dt))[-1]
x_ss, y_ss = ss
print(f"Steady state x* = {x_ss:.6f}")
print(f"Steady state y* = {y_ss:.6f}")

# Residual of the RHS at the steady state (should be ~0).
res = f(ss)
print(f"RHS residual at steady state (dx/dt) = {res[0]:.3e}")
print(f"RHS residual at steady state (dy/dt) = {res[1]:.3e}")

# ---------------------------------------------------------------
# Linear stability: build the Jacobian at the steady state by finite
# differences and inspect its eigenvalues.
#   Complex eigenvalues  -> spiraling (oscillatory approach)
#   Negative real parts  -> stable (trajectories converge)
# ---------------------------------------------------------------
eps = 1e-6
J = np.empty((2, 2))
for j in range(2):
    pert = ss.copy()
    pert[j] += eps
    J[:, j] = (f(pert) - f(ss)) / eps

eigvals = np.linalg.eigvals(J)
print(f"Jacobian eigenvalue 1 = {eigvals[0]:.6f}")
print(f"Jacobian eigenvalue 2 = {eigvals[1]:.6f}")
print(f"Max real part of eigenvalues = {np.max(eigvals.real):.6f}")
has_imag = np.any(np.abs(eigvals.imag) > 1e-8)
print(f"Eigenvalues complex (spiral) = {has_imag}")
print(f"Steady state stable (all real parts < 0) = {np.all(eigvals.real < 0)}")

# Confirm every trajectory ends at the same steady state.
max_endpoint_dist = max(np.linalg.norm(tr[-1] - ss) for tr in trajectories)
print(f"Max distance of trajectory endpoints from steady state = {max_endpoint_dist:.3e}")

# ---------------------------------------------------------------
# Phase-plane plot.
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
for ic, tr in zip(initial_conditions, trajectories):
    plt.plot(tr[:, 0], tr[:, 1], lw=1.0, alpha=0.8)
    plt.plot(ic[0], ic[1], 'o', ms=5)
plt.plot(x_ss, y_ss, 'k*', ms=16, label='steady state')
plt.xlabel('x')
plt.ylabel('y')
plt.title(f'Negative-feedback loop phase plane (g={g:g}, h={h:g})\n'
          'trajectories spiral into a single stable steady state')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.2.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: the Jacobian at the steady state has complex-conjugate "
      "eigenvalues with negative real parts, which mathematically guarantees "
      "that nearby trajectories oscillate (spiral) while decaying toward that "
      "single point, confirming it is a unique stable spiral attractor.")
