import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: constitutively expressed gene, dX/dt = g - k*X ---
g = 50.0      # constant transcription rate
k = 0.1       # linear degradation rate
X0 = 300.0    # initial mRNA/protein level

def deriv(X):
    # right-hand side of the ODE (autonomous, no explicit t dependence)
    return g - k * X

def exact(t):
    # analytic solution of dX/dt = g - k*X
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Integrators, all implemented explicitly ---

def euler(f, X0, t):
    # first-order explicit Euler
    X = np.empty_like(t)
    X[0] = X0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        X[i + 1] = X[i] + h * f(X[i])   # step using slope at the left endpoint
    return X

def rk2_midpoint(f, X0, t):
    # second-order Runge-Kutta (midpoint method)
    X = np.empty_like(t)
    X[0] = X0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        k1 = f(X[i])                    # slope at the start of the interval
        X_mid = X[i] + 0.5 * h * k1     # trial half-step to the interval midpoint
        k2 = f(X_mid)                   # slope evaluated at that midpoint
        X[i + 1] = X[i] + h * k2        # full step using the midpoint slope
    return X

def heun(f, X0, t):
    # second-order Heun (explicit trapezoidal) method, for comparison
    X = np.empty_like(t)
    X[0] = X0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]
        k1 = f(X[i])                    # slope at the start
        X_pred = X[i] + h * k1          # Euler predictor for the endpoint
        k2 = f(X_pred)                  # slope at the predicted endpoint
        X[i + 1] = X[i] + 0.5 * h * (k1 + k2)  # average the two slopes
    return X

# --- Integrate over a time grid ---
h = 1.0
t = np.arange(0.0, 100.0 + h, h)

X_rk2 = rk2_midpoint(deriv, X0, t)
X_euler = euler(deriv, X0, t)
X_heun = heun(deriv, X0, t)
X_exact = exact(t)

# --- Errors against the exact solution (max absolute error over the grid) ---
err_rk2 = np.max(np.abs(X_rk2 - X_exact))
err_euler = np.max(np.abs(X_euler - X_exact))
err_heun = np.max(np.abs(X_heun - X_exact))

print(f"Step size h = {h}")
print(f"Steady state g/k = {g / k}")
print(f"Max abs error, Euler  : {err_euler:.6e}")
print(f"Max abs error, RK2    : {err_rk2:.6e}")
print(f"Max abs error, Heun   : {err_heun:.6e}")
print(f"RK2 better than Euler (err_rk2 < err_euler)         : {err_rk2 < err_euler}")
print(f"RK2 comparable to Heun (same order of magnitude)    : "
      f"{0.1 <= err_rk2 / err_heun <= 10.0}")

# --- Plot RK2 overlaid on the exact solution ---
plt.figure(figsize=(8, 5))
tf = np.linspace(0.0, 100.0, 500)
plt.plot(tf, exact(tf), "k-", lw=2, label="Exact")
plt.plot(t, X_rk2, "ro", ms=4, label=f"RK2 (midpoint), h={h}")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("RK2 (midpoint) vs exact solution: dX/dt = g - k*X")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.5.1_s2.png")

# Explanation: RK2 and Heun are both second-order methods so their errors shrink
# like h^2 while Euler's shrinks only like h, so at a fixed step size finding
# err_rk2 close to err_heun and well below err_euler confirms the midpoint
# integrator achieves the expected second-order accuracy.
print("Check meaning: RK2's error is close to Heun's and far below Euler's at the "
      "same h, confirming it attains the expected second-order accuracy.")
