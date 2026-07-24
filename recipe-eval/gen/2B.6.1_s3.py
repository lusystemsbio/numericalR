import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model: constitutive gene expression dX/dt = g - k*X -----
g, k, X0 = 50.0, 0.1, 300.0

def f(t, X):
    # RHS of the ODE: production g minus linear degradation k*X
    return g - k * X

def exact(t):
    # Analytic solution of the linear ODE
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ----- Explicit fourth-order Runge-Kutta step -----
def rk4_step(t, X, h):
    k1 = f(t, X)                    # slope at the start of the step
    k2 = f(t + 0.5 * h, X + 0.5 * h * k1)  # slope at the midpoint (using k1)
    k3 = f(t + 0.5 * h, X + 0.5 * h * k2)  # slope again at midpoint (using k2)
    k4 = f(t + h, X + h * k3)       # slope at the end of the step (using k3)
    # weighted average of the four slopes: weights 1, 2, 2, 1 (sum 6)
    return X + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)

# ----- Second-order Runge-Kutta (midpoint) step, for comparison -----
def rk2_step(t, X, h):
    k1 = f(t, X)                    # slope at the start
    k2 = f(t + 0.5 * h, X + 0.5 * h * k1)  # slope at the midpoint
    return X + h * k2               # advance using the midpoint slope

# ----- Explicit (forward) Euler step, for comparison -----
def euler_step(t, X, h):
    return X + h * f(t, X)          # advance using only the start slope

# ----- Generic fixed-step integrator driver -----
def integrate(step, h, t_end):
    n = int(round(t_end / h))
    t = np.linspace(0.0, n * h, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        X[i + 1] = step(t[i], X[i], h)
    return t, X

# ----- Run RK4 and generate the overlay plot -----
t_end = 50.0
h_plot = 1.0
t_rk4, X_rk4 = integrate(rk4_step, h_plot, t_end)
t_fine = np.linspace(0.0, t_end, 1000)

plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine), 'k-', lw=2, label="Exact")
plt.plot(t_rk4, X_rk4, 'ro', ms=5, label="RK4 (h=%.1f)" % h_plot)
plt.xlabel("time")
plt.ylabel("X(t)")
plt.title("RK4 vs exact: dX/dt = g - k*X  (g=%.0f, k=%.1f, X0=%.0f)" % (g, k, X0))
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.6.1_s3.png")

# ----- Accuracy check: global error at t_end for each method -----
def max_error(step, h):
    t, X = integrate(step, h, t_end)
    return np.max(np.abs(X - exact(t)))

h1 = 1.0
h2 = 0.5

err_euler = max_error(euler_step, h1)
err_rk2 = max_error(rk2_step, h1)
err_rk4_h1 = max_error(rk4_step, h1)
err_rk4_h2 = max_error(rk4_step, h2)

print("Max error, Euler (h=%.2f):        %.6e" % (h1, err_euler))
print("Max error, RK2   (h=%.2f):        %.6e" % (h1, err_rk2))
print("Max error, RK4   (h=%.2f):        %.6e" % (h1, err_rk4_h1))
print("Max error, RK4   (h=%.2f):        %.6e" % (h2, err_rk4_h2))
print("RK4 more accurate than RK2 by factor: %.2f" % (err_rk2 / err_rk4_h1))
print("RK4 more accurate than Euler by factor: %.2f" % (err_euler / err_rk4_h1))
print("RK4 error ratio when step halved (expect ~16): %.2f" % (err_rk4_h1 / err_rk4_h2))

# Explanation:
# A 4th-order method has global error O(h^4), so halving h should reduce the
# error by 2^4 = 16; observing that ~16x drop (and errors far below RK2/Euler)
# confirms the integrator truly achieves fourth-order accuracy.
print("Explanation: RK4 has O(h^4) global error, so halving the step should cut "
      "the error by 2^4=16; seeing that ~16x reduction confirms fourth-order accuracy.")
