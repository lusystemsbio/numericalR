import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g = 50.0    # constant transcription rate
k = 0.1     # linear degradation rate
X0 = 300.0  # initial mRNA/protein level

# --- ODE right-hand side: dX/dt = g - k*X ---
def f(X):
    return g - k * X

# --- Exact analytic solution ---
def exact(t):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Second-order Runge-Kutta (midpoint) integrator, done explicitly ---
def rK2(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + h, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        X = Xs[i]
        k1 = f(X)                    # slope at the start of the interval
        X_mid = X + 0.5 * h * k1     # trial half-step to the interval midpoint
        k2 = f(X_mid)                # slope evaluated at that midpoint
        Xs[i + 1] = X + h * k2       # full step using the midpoint slope
    return ts, Xs

# --- Explicit Heun (RK2 trapezoidal) for comparison ---
def heun(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + h, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        X = Xs[i]
        k1 = f(X)                    # slope at the start
        X_pred = X + h * k1          # Euler predictor for the endpoint
        k2 = f(X_pred)               # slope at the predicted endpoint
        Xs[i + 1] = X + 0.5 * h * (k1 + k2)  # average the two slopes
    return ts, Xs

# --- Explicit forward Euler for comparison ---
def euler(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + h, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        Xs[i + 1] = Xs[i] + h * f(Xs[i])  # single slope at the start
    return ts, Xs

# --- Integrate over a fixed interval with a common step size ---
t0, tf, h = 0.0, 60.0, 1.0
t_rk2, X_rk2 = rK2(f, X0, t0, tf, h)
t_heun, X_heun = heun(f, X0, t0, tf, h)
t_eul, X_eul = euler(f, X0, t0, tf, h)
X_ex = exact(t_rk2)

# --- Max absolute error against the exact solution ---
err_rk2 = np.max(np.abs(X_rk2 - X_ex))
err_heun = np.max(np.abs(X_heun - exact(t_heun)))
err_eul = np.max(np.abs(X_eul - exact(t_eul)))

print(f"Step size h                = {h}")
print(f"Max abs error RK2 (midpoint) = {err_rk2:.6e}")
print(f"Max abs error Heun           = {err_heun:.6e}")
print(f"Max abs error Euler          = {err_eul:.6e}")
print(f"RK2 error / Heun error       = {err_rk2 / err_heun:.6f}")
print(f"RK2 error / Euler error      = {err_rk2 / err_eul:.6f}")
print(f"Euler error / RK2 error      = {err_eul / err_rk2:.6f}")

# Explanation: RK2's error being of the same order as Heun's yet much smaller
# than Euler's at the identical step size confirms RK2 is a genuine second-order
# method, since a true O(h^2) scheme must beat the O(h) Euler and match the other
# second-order method when compared to the exact solution.
print("Explanation: RK2 matches Heun (both O(h^2)) and is far below Euler (O(h)) "
      "at the same h, which is the signature of a correct second-order method.")

# --- Plot RK2 overlaid on the exact solution ---
t_fine = np.linspace(t0, tf, 500)
plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine), 'k-', label='Exact')
plt.plot(t_rk2, X_rk2, 'ro', markersize=4, label='RK2 (midpoint)')
plt.axhline(g / k, color='gray', ls='--', lw=0.8, label='Steady state g/k')
plt.xlabel('time t')
plt.ylabel('X(t)')
plt.title('Constitutive gene expression: RK2 (midpoint) vs exact')
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.5.1_s3.png")
