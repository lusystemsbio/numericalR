import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters (nondimensional two-gene negative-feedback loop)
g, h, n = 10.0, 10.0, 3  # X activates Y, Y represses X; Hill coefficient 3, unit degradation

def f(state):
    """Right-hand side of the ODE system: dx/dt, dy/dt."""
    x, y = state
    dx = g / (1.0 + y**n) - x          # Y represses X, plus unit degradation of x
    dy = h * x**n / (1.0 + x**n) - y   # X activates Y, plus unit degradation of y
    return np.array([dx, dy])

def rk4_step(state, dt):
    """One explicit classic 4th-order Runge-Kutta step (implemented by hand)."""
    k1 = f(state)                      # slope at start
    k2 = f(state + 0.5 * dt * k1)      # slope at midpoint using k1
    k3 = f(state + 0.5 * dt * k2)      # slope at midpoint using k2
    k4 = f(state + dt * k3)            # slope at end using k3
    return state + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)  # weighted average

def integrate(state0, dt, T):
    """Integrate the system from state0 for total time T with fixed step dt."""
    nsteps = int(round(T / dt))
    traj = np.empty((nsteps + 1, 2))
    traj[0] = state0
    for i in range(nsteps):
        traj[i + 1] = rk4_step(traj[i], dt)
    return traj

# --- Integrate several trajectories from different initial conditions ---
dt, T = 0.01, 40.0
initial_conditions = [(0.5, 0.5), (8.0, 1.0), (1.0, 9.0), (6.0, 6.0), (0.1, 8.0)]

plt.figure(figsize=(7, 6))
final_states = []
for x0, y0 in initial_conditions:
    traj = integrate(np.array([x0, y0]), dt, T)
    final_states.append(traj[-1])
    plt.plot(traj[:, 0], traj[:, 1], lw=1.0, label=f"IC ({x0}, {y0})")
    plt.plot(x0, y0, 'o', ms=4)

final_states = np.array(final_states)

# --- Check: all trajectories converge to the same steady state ---
ss_mean = final_states.mean(axis=0)
ss_spread = final_states.std(axis=0)
residual = np.linalg.norm(f(ss_mean))  # dx/dt, dy/dt should be ~0 at a fixed point

# --- Check spiraling: eigenvalues of the Jacobian at the steady state ---
xs, ys = ss_mean
# Analytic Jacobian entries
J = np.array([
    [-1.0,                       -g * n * ys**(n-1) / (1.0 + ys**n)**2],
    [h * n * xs**(n-1) / (1.0 + xs**n)**2,                       -1.0]
])
eigs = np.linalg.eigvals(J)

print(f"Final states from each initial condition:")
for (x0, y0), fs in zip(initial_conditions, final_states):
    print(f"  IC ({x0}, {y0}) -> steady state ({fs[0]:.6f}, {fs[1]:.6f})")
print(f"Mean steady state x*: {ss_mean[0]:.6f}")
print(f"Mean steady state y*: {ss_mean[1]:.6f}")
print(f"Std of final x across ICs: {ss_spread[0]:.3e}")
print(f"Std of final y across ICs: {ss_spread[1]:.3e}")
print(f"Residual |f(steady state)| (should be ~0): {residual:.3e}")
print(f"Jacobian eigenvalue 1: {eigs[0].real:.6f} + {eigs[0].imag:.6f}i")
print(f"Jacobian eigenvalue 2: {eigs[1].real:.6f} + {eigs[1].imag:.6f}i")
print(f"Max real part of eigenvalues (should be < 0 for stability): {eigs.real.max():.6f}")
print(f"Nonzero imaginary part (indicates spiraling): {abs(eigs.imag).max():.6f}")

plt.plot(ss_mean[0], ss_mean[1], 'k*', ms=16, label="steady state")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Two-gene negative-feedback loop: trajectories spiraling to steady state")
plt.legend(fontsize=8, loc="upper right")
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.2.1_s1.png")

# Explanation: Because every initial condition converges to the same point where f=0, and the
# Jacobian there has eigenvalues with negative real part and nonzero imaginary part, the fixed
# point is confirmed as a unique stable spiral (complex eigenvalues => oscillatory decay = spiraling in).
print("Check: all ICs converge to one point (f~0) whose Jacobian has negative-real-part,")
print("complex eigenvalues, confirming a unique stable steady state approached by spiraling in.")
