import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Model definition: dX/dt = g - k*X ---
g = 50.0
k = 0.1
X0 = 300.0

def f(X):
    # RHS of the ODE (right-hand side / slope function)
    return g - k * X

def exact(t):
    # Analytic solution X(t) = g/k + (X0 - g/k)*exp(-k*t)
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Explicit Heun (2nd-order) integrator ---
def heun(h, t_end):
    n = int(round(t_end / h))
    t = np.linspace(0.0, n * h, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        s1 = f(X[i])                 # slope at the start of the step
        X_pred = X[i] + h * s1       # Euler predictor: step to a tentative endpoint
        s2 = f(X_pred)               # slope evaluated at the predicted endpoint
        X[i + 1] = X[i] + h * 0.5 * (s1 + s2)  # advance by the average of the two slopes
    return t, X

# --- Explicit forward Euler, for comparison ---
def euler(h, t_end):
    n = int(round(t_end / h))
    t = np.linspace(0.0, n * h, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        X[i + 1] = X[i] + h * f(X[i])  # single slope, evaluated only at the start
    return t, X

t_end = 60.0

# --- Solve and plot ---
h = 2.0
t_h, X_h = heun(h, t_end)
t_e, X_e = euler(h, t_end)
X_exact_on_grid = exact(t_h)

t_fine = np.linspace(0.0, t_end, 500)
plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine), "k-", label="Exact")
plt.plot(t_h, X_h, "o--", color="tab:blue", label=f"Heun (h={h})")
plt.plot(t_e, X_e, "s--", color="tab:red", alpha=0.6, label=f"Euler (h={h})")
plt.xlabel("t")
plt.ylabel("X(t)")
plt.title("Constitutive gene expression: Heun vs exact")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.4.1_s3.png")

# --- Check 1: Heun beats Euler at the same step size ---
# Use max absolute error over the shared grid.
err_heun = np.max(np.abs(X_h - X_exact_on_grid))
err_euler = np.max(np.abs(X_e - exact(t_e)))
print(f"Max error, Euler (h={h}): {err_euler:.6e}")
print(f"Max error, Heun  (h={h}): {err_heun:.6e}")
print(f"Heun more accurate than Euler at same h: {err_heun < err_euler}")

# --- Check 2: halving h cuts Heun's error by ~4x (2nd order: error ~ h^2) ---
h2 = h / 2.0
t_h2, X_h2 = heun(h2, t_end)
err_heun_h2 = np.max(np.abs(X_h2 - exact(t_h2)))
ratio = err_heun / err_heun_h2
print(f"Max error, Heun  (h={h2}): {err_heun_h2:.6e}")
print(f"Error ratio (h vs h/2): {ratio:.4f}")
print(f"Observed order of convergence (log2 of ratio): {np.log2(ratio):.4f}")

# One-sentence explanation:
# The error dropping by ~4x when h is halved means error scales like h^2,
# which is exactly the signature of a second-order method, confirming Heun
# is genuinely second-order accurate (as opposed to Euler's first-order, ~2x drop).
print("Explanation: an error ratio near 4 (order ~2) confirms the scheme is "
      "second-order accurate, since error ~ h^2 means halving h reduces error fourfold.")
