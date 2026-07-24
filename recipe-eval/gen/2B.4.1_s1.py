import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# ODE right-hand side: dX/dt = g - k*X
def rhs(X, g, k):
    return g - k * X


# Exact analytic solution
def exact(t, g, k, X0):
    return g / k + (X0 - g / k) * np.exp(-k * t)


# Explicit second-order Heun integrator
def heun(g, k, X0, t_end, dt):
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        f0 = rhs(X[i], g, k)              # slope at the start of the step
        X_pred = X[i] + dt * f0           # Euler predictor to the endpoint
        f1 = rhs(X_pred, g, k)            # slope at the predicted endpoint
        X[i + 1] = X[i] + dt * 0.5 * (f0 + f1)  # advance by the average slope
    return t, X


# Explicit forward Euler integrator (for comparison)
def euler(g, k, X0, t_end, dt):
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        X[i + 1] = X[i] + dt * rhs(X[i], g, k)  # single slope at the start
    return t, X


# Parameters
g, k, X0 = 50.0, 0.1, 300.0
t_end = 50.0

# Max-absolute-error helper against the exact solution
def max_err(t, X):
    return np.max(np.abs(X - exact(t, g, k, X0)))

# Reference step size
dt = 1.0
t_h, X_h = heun(g, k, X0, t_end, dt)
t_e, X_e = euler(g, k, X0, t_end, dt)

err_heun = max_err(t_h, X_h)
err_euler = max_err(t_e, X_e)

# Halved step size for Heun to check order of accuracy
dt2 = 0.5
t_h2, X_h2 = heun(g, k, X0, t_end, dt2)
err_heun_half = max_err(t_h2, X_h2)
ratio = err_heun / err_heun_half

# Reported numerical results
print(f"Step size dt                       = {dt}")
print(f"Heun  max abs error (dt)           = {err_heun:.6e}")
print(f"Euler max abs error (dt)           = {err_euler:.6e}")
print(f"Heun is more accurate than Euler   = {err_heun < err_euler}")
print(f"Heun max abs error (dt/2 = {dt2})    = {err_heun_half:.6e}")
print(f"Error reduction ratio (dt / dt/2)  = {ratio:.4f}")
print("Expected ratio for 2nd-order method ~ 4")

# Plot: Heun overlaid on the exact solution
t_fine = np.linspace(0.0, t_end, 500)
plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine, g, k, X0), 'k-', lw=2, label="Exact")
plt.plot(t_h, X_h, 'ro--', ms=4, label=f"Heun (dt={dt})")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Constitutive gene expression: Heun vs exact")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.4.1_s1.png")

# Explanation:
# The check confirms the result because a lower error than Euler at the same
# step, together with the ~4x error drop when the step is halved, is exactly
# the O(dt^2) local behavior a correct second-order Heun method must exhibit.
print("Explanation: beating Euler at equal dt and cutting error ~4x when dt "
      "is halved is the O(dt^2) signature that confirms a correct 2nd-order Heun method.")
