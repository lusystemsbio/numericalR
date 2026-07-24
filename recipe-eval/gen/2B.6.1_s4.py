import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ODE right-hand side: dX/dt = g - k*X
def f(t, X, g, k):
    return g - k * X


# Exact analytical solution for comparison
def exact(t, g, k, X0):
    return g / k + (X0 - g / k) * np.exp(-k * t)


# --- Explicit fourth-order Runge-Kutta integrator ---
def rk4(f, X0, t0, tf, h, g, k):
    n = int(round((tf - t0) / h))          # number of steps
    ts = np.empty(n + 1)
    Xs = np.empty(n + 1)
    ts[0], Xs[0] = t0, X0
    t, X = t0, X0
    for i in range(n):
        k1 = f(t, X, g, k)                  # slope at the start of the step
        k2 = f(t + h / 2, X + h / 2 * k1, g, k)  # slope at the midpoint using k1
        k3 = f(t + h / 2, X + h / 2 * k2, g, k)  # slope at the midpoint using k2
        k4 = f(t + h, X + h * k3, g, k)     # slope at the end using k3
        # combine with weights 1, 2, 2, 1 (divide by 6)
        X = X + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
        ts[i + 1], Xs[i + 1] = t, X
    return ts, Xs


# --- Euler and RK2 (midpoint) for the accuracy comparison ---
def euler(f, X0, t0, tf, h, g, k):
    n = int(round((tf - t0) / h))
    ts = np.empty(n + 1); Xs = np.empty(n + 1)
    ts[0], Xs[0] = t0, X0
    t, X = t0, X0
    for i in range(n):
        X = X + h * f(t, X, g, k)           # single forward slope
        t = t + h
        ts[i + 1], Xs[i + 1] = t, X
    return ts, Xs


def rk2(f, X0, t0, tf, h, g, k):
    n = int(round((tf - t0) / h))
    ts = np.empty(n + 1); Xs = np.empty(n + 1)
    ts[0], Xs[0] = t0, X0
    t, X = t0, X0
    for i in range(n):
        k1 = f(t, X, g, k)                  # start slope
        k2 = f(t + h / 2, X + h / 2 * k1, g, k)  # midpoint slope
        X = X + h * k2                      # step with midpoint slope
        t = t + h
        ts[i + 1], Xs[i + 1] = t, X
    return ts, Xs


# Parameters
g, k, X0 = 50.0, 0.1, 300.0
t0, tf = 0.0, 100.0
h = 1.0

# Integrate with RK4 and evaluate exact solution on the same grid
t_rk4, X_rk4 = rk4(f, X0, t0, tf, h, g, k)
X_exact = exact(t_rk4, g, k, X0)

# Steady-state value g/k for reference
print(f"Steady-state value g/k: {g / k:.6f}")
print(f"Final RK4 value at t={tf}: {X_rk4[-1]:.10f}")
print(f"Final exact value at t={tf}: {X_exact[-1]:.10f}")

# Max absolute error of each method at step size h
_, Xe = euler(f, X0, t0, tf, h, g, k)
_, X2 = rk2(f, X0, t0, tf, h, g, k)
te = exact(t_rk4, g, k, X0)
err_euler = np.max(np.abs(Xe - te))
err_rk2 = np.max(np.abs(X2 - te))
err_rk4 = np.max(np.abs(X_rk4 - te))
print(f"Max error Euler (h={h}): {err_euler:.6e}")
print(f"Max error RK2   (h={h}): {err_rk2:.6e}")
print(f"Max error RK4   (h={h}): {err_rk4:.6e}")

# Convergence check: halve the step, error should drop ~16x (4th order)
h2 = h / 2
t_rk4_h2, X_rk4_h2 = rk4(f, X0, t0, tf, h2, g, k)
err_rk4_h2 = np.max(np.abs(X_rk4_h2 - exact(t_rk4_h2, g, k, X0)))
print(f"Max error RK4   (h={h2}): {err_rk4_h2:.6e}")
ratio = err_rk4 / err_rk4_h2
print(f"Error ratio (h -> h/2): {ratio:.4f}  (expected ~16 for 4th order)")

# Plot RK4 overlaid on exact solution
plt.figure(figsize=(9, 5))
plt.plot(t_rk4, X_exact, 'b-', lw=2, label='Exact')
plt.plot(t_rk4, X_rk4, 'ro', ms=4, label=f'RK4 (h={h})')
plt.axhline(g / k, color='gray', ls='--', lw=1, label='Steady state g/k')
plt.xlabel('time t')
plt.ylabel('X(t)')
plt.title('Constitutive gene expression: RK4 vs exact')
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.6.1_s4.png")

# One-sentence explanation of why the convergence check confirms correctness:
print("Explanation: An error ratio near 16 when the step is halved means the "
      "error scales as h^4, which is the defining signature of a correctly "
      "implemented fourth-order method (2^4 = 16), so RK4's far smaller error "
      "than Euler/RK2 plus this scaling confirms the integrator is right.")
