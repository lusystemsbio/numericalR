import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Model: constitutively expressed gene, dX/dt = g - k*X ---
g, k, X0 = 50.0, 0.1, 300.0

def f(X):
    # right-hand side of the ODE (autonomous here)
    return g - k * X

def exact(t):
    # analytic solution
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Integrators (all explicit, one step per call) ---
def step_euler(X, h):
    # forward Euler: advance using slope at the start of the interval
    return X + h * f(X)

def step_heun(X, h):
    k1 = f(X)                 # slope at start
    X_pred = X + h * k1        # Euler predictor to end of interval
    k2 = f(X_pred)            # slope at predicted end
    return X + h * 0.5 * (k1 + k2)   # average the two slopes

def step_rk2_midpoint(X, h):
    k1 = f(X)                     # slope at the start of the interval
    X_mid = X + 0.5 * h * k1       # trial half-step to the midpoint
    k2 = f(X_mid)                 # slope evaluated at the midpoint
    return X + h * k2              # full step using the midpoint slope

def integrate(step, t_end, h):
    n = int(round(t_end / h))
    t = np.linspace(0.0, n * h, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        X[i + 1] = step(X[i], h)
    return t, X

# --- Run RK2 (midpoint) and compare against exact solution ---
t_end, h = 30.0, 1.0
t, X_rk2 = integrate(step_rk2_midpoint, t_end, h)
X_ex = exact(t)

# --- Plot: RK2 overlaid on exact ---
tt = np.linspace(0.0, t_end, 400)
plt.figure(figsize=(8, 5))
plt.plot(tt, exact(tt), "k-", label="Exact")
plt.plot(t, X_rk2, "ro", ms=5, label="RK2 (midpoint), h=%.2g" % h)
plt.xlabel("time")
plt.ylabel("X(t)")
plt.title("RK2 (midpoint) vs exact:  dX/dt = g - k*X")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.5.1_s1.png")

# --- Accuracy check at the same step size ---
_, X_eu = integrate(step_euler, t_end, h)
_, X_he = integrate(step_heun, t_end, h)

def max_abs_err(X):
    return np.max(np.abs(X - X_ex))

err_euler = max_abs_err(X_eu)
err_heun  = max_abs_err(X_he)
err_rk2   = max_abs_err(X_rk2)

print("Parameters: g = %.4g, k = %.4g, X0 = %.4g, h = %.4g" % (g, k, X0, h))
print("X(t_end) exact           = %.6f" % X_ex[-1])
print("X(t_end) RK2 (midpoint)  = %.6f" % X_rk2[-1])
print("Max abs error, Euler     = %.6e" % err_euler)
print("Max abs error, Heun      = %.6e" % err_heun)
print("Max abs error, RK2       = %.6e" % err_rk2)
print("RK2 vs Heun error ratio  = %.6f" % (err_rk2 / err_heun))
print("Euler / RK2 error ratio  = %.6f" % (err_euler / err_rk2))
print("RK2 comparable to Heun (within 2x)?  %s" % (0.5 < err_rk2 / err_heun < 2.0))
print("RK2 better than Euler?               %s" % (err_rk2 < err_euler))

# One-sentence explanation:
print("Explanation: Because both RK2 variants are second-order (global error ~ h^2) "
      "while Euler is only first-order (~ h), finding the RK2 error on par with Heun and "
      "well below Euler's at the same h confirms the midpoint step was implemented with the "
      "correct second-order accuracy.")
