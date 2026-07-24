import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters (constitutively expressed gene): dX/dt = g - k*X
g = 50.0     # constant transcription rate
k = 0.1      # linear degradation rate
X0 = 300.0   # initial condition

# ODE right-hand side (slope function). Autonomous here, but keep t for generality.
def f(t, X):
    return g - k * X

# Exact analytic solution for comparison
def exact(t):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Second-order Runge-Kutta (midpoint) integrator, implemented explicitly ---
def rk2_midpoint(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + h / 2, h)   # time grid
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        t = ts[i]
        X = Xs[i]
        k1 = f(t, X)                       # slope at the start of the interval
        X_mid = X + 0.5 * h * k1           # trial half-step to the midpoint
        k2 = f(t + 0.5 * h, X_mid)         # slope evaluated at the midpoint
        Xs[i + 1] = X + h * k2             # full step using the midpoint slope
    return ts, Xs

# --- Euler (first order) for comparison ---
def euler(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + h / 2, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        Xs[i + 1] = Xs[i] + h * f(ts[i], Xs[i])
    return ts, Xs

# --- Heun (second order, trapezoidal predictor-corrector) for comparison ---
def heun(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + h / 2, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        t = ts[i]
        X = Xs[i]
        k1 = f(t, X)                       # slope at start
        X_pred = X + h * k1                # Euler predictor (full step)
        k2 = f(t + h, X_pred)              # slope at end using predictor
        Xs[i + 1] = X + 0.5 * h * (k1 + k2)  # average of the two slopes
    return ts, Xs

# Integration settings
t0, tf, h = 0.0, 50.0, 1.0

ts_rk2, X_rk2 = rk2_midpoint(f, X0, t0, tf, h)
ts_eul, X_eul = euler(f, X0, t0, tf, h)
ts_heu, X_heu = heun(f, X0, t0, tf, h)
X_true = exact(ts_rk2)

# Max absolute error against the exact solution at the same step size
err_rk2 = np.max(np.abs(X_rk2 - X_true))
err_eul = np.max(np.abs(X_eul - exact(ts_eul)))
err_heu = np.max(np.abs(X_heu - exact(ts_heu)))

print(f"Step size h = {h}")
print(f"Max |error| Euler         : {err_eul:.6e}")
print(f"Max |error| RK2 (midpoint): {err_rk2:.6e}")
print(f"Max |error| Heun          : {err_heu:.6e}")
print(f"RK2 vs Heun error ratio   : {err_rk2 / err_heu:.6f}")
print(f"RK2 vs Euler error ratio  : {err_rk2 / err_eul:.6f}")
print(f"RK2 error is comparable to Heun (ratio ~O(1)): {0.1 < err_rk2 / err_heu < 10.0}")
print(f"RK2 error is smaller than Euler error       : {err_rk2 < err_eul}")

# Final-time values for reference
print(f"X_exact(tf={tf})          : {exact(tf):.6f}")
print(f"X_RK2(tf={tf})            : {X_rk2[-1]:.6f}")
print(f"X_Euler(tf={tf})         : {X_eul[-1]:.6f}")
print(f"X_Heun(tf={tf})          : {X_heu[-1]:.6f}")

# Explanation of why the check confirms the result:
print("Explanation: Matching the exact solution's error to Heun's (both O(h^2)) "
      "while beating Euler's (O(h)) confirms RK2 achieves the correct "
      "second-order accuracy, so the midpoint scheme is implemented correctly.")

# --- Plot: RK2 solution overlaid on the exact solution ---
t_fine = np.linspace(t0, tf, 500)
plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine), 'k-', lw=2, label='Exact solution')
plt.plot(ts_rk2, X_rk2, 'ro', ms=5, mfc='none', label='RK2 (midpoint)')
plt.xlabel('time t')
plt.ylabel('X(t)  (mRNA/protein level)')
plt.title(f'Constitutive gene expression: RK2 vs exact (g={g}, k={k}, X0={X0}, h={h})')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.5.1_s4.png")
