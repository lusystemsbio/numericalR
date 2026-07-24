import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model definition: constitutively expressed gene ----
# dX/dt = g - k*X
g = 50.0     # constant transcription rate
k = 0.1      # linear degradation rate
X0 = 300.0   # initial amount

def f(X):
    # right-hand side of the ODE
    return g - k * X

def exact(t):
    # analytic solution X(t) = g/k + (X0 - g/k)*exp(-k*t)
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ---- Explicit Euler integrator (for comparison) ----
def euler(h, T):
    n = int(round(T / h))
    ts = np.linspace(0.0, n * h, n + 1)
    Xs = np.empty(n + 1)
    Xs[0] = X0
    for i in range(n):
        Xs[i + 1] = Xs[i] + h * f(Xs[i])   # advance by the starting slope
    return ts, Xs

# ---- Explicit Heun (second-order predictor-corrector) integrator ----
def heun(h, T):
    n = int(round(T / h))
    ts = np.linspace(0.0, n * h, n + 1)
    Xs = np.empty(n + 1)
    Xs[0] = X0
    for i in range(n):
        s1 = f(Xs[i])                       # slope at the start of the step
        X_pred = Xs[i] + h * s1             # Euler step to a predicted endpoint
        s2 = f(X_pred)                      # slope evaluated at the predicted endpoint
        Xs[i + 1] = Xs[i] + h * 0.5 * (s1 + s2)  # advance by the average of the two slopes
    return ts, Xs

# ---- Run the integrators ----
T = 60.0     # total time
h = 2.0      # step size

t_heun, X_heun = heun(h, T)
t_euler, X_euler = euler(h, T)

# Dense exact curve for plotting, and exact values at the integrator grid points
t_dense = np.linspace(0.0, T, 500)
X_dense = exact(t_dense)
X_exact_grid = exact(t_heun)

# ---- Plot: Heun overlaid on the exact solution ----
plt.figure(figsize=(8, 5))
plt.plot(t_dense, X_dense, "k-", lw=2, label="Exact")
plt.plot(t_heun, X_heun, "ro--", ms=5, label=f"Heun (h={h})")
plt.plot(t_euler, X_euler, "b^:", ms=5, label=f"Euler (h={h})")
plt.axhline(g / k, color="gray", ls=":", lw=1, label="Steady state g/k")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Constitutive gene expression: Heun vs Euler vs Exact")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.4.1_s5.png")

# ---- Check 1: at the same step size, Heun beats Euler ----
# Use the maximum absolute error over all grid points.
err_heun = np.max(np.abs(X_heun - X_exact_grid))
err_euler = np.max(np.abs(X_euler - X_exact_grid))

print(f"Step size h                          : {h}")
print(f"Max abs error, Euler (h={h})         : {err_euler:.6e}")
print(f"Max abs error, Heun  (h={h})         : {err_heun:.6e}")
print(f"Heun is closer than Euler            : {err_heun < err_euler}")

# ---- Check 2: halving the step should cut Heun's error ~4x (second order) ----
_, X_heun_h  = heun(h, T)
_, X_heun_h2 = heun(h / 2.0, T)

err_h  = np.max(np.abs(X_heun_h  - exact(np.linspace(0.0, T, int(round(T / h)) + 1))))
err_h2 = np.max(np.abs(X_heun_h2 - exact(np.linspace(0.0, T, int(round(T / (h / 2.0))) + 1))))
ratio = err_h / err_h2

print(f"Max abs error, Heun (h={h})          : {err_h:.6e}")
print(f"Max abs error, Heun (h={h/2.0})        : {err_h2:.6e}")
print(f"Error reduction ratio (should ~4)    : {ratio:.4f}")

# Explanation: A second-order method has global error proportional to h^2, so
# halving h divides the error by 2^2 = 4 -- observing the ~4x drop (and Heun
# beating Euler at equal h) confirms the integrator is genuinely second order.
print("Explanation: global error ~ h^2 for a 2nd-order method, so halving h "
      "gives a 2^2=4x error drop, confirming Heun's second-order accuracy.")
