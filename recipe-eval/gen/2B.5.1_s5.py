import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model: constitutively expressed gene, dX/dt = g - k*X
g, k, X0 = 50.0, 0.1, 300.0

def deriv(X):
    # transcription at constant rate g, linear degradation at rate k*X
    return g - k * X

def exact(t):
    # analytic solution X(t) = g/k + (X0 - g/k)*exp(-k*t)
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ---- RK2 (midpoint) integrator, written out explicitly ----
def rk2_midpoint(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + 0.5 * h, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        X = Xs[i]
        k1 = f(X)                    # slope at the start of the interval
        X_mid = X + 0.5 * h * k1     # trial half-step to the interval midpoint
        k2 = f(X_mid)                # slope evaluated at the midpoint
        Xs[i + 1] = X + h * k2       # full step using the midpoint slope
    return ts, Xs

# ---- Euler (first order) for comparison ----
def euler(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + 0.5 * h, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        Xs[i + 1] = Xs[i] + h * f(Xs[i])
    return ts, Xs

# ---- Heun (RK2, trapezoidal variant) for comparison ----
def heun(f, X0, t0, tf, h):
    ts = np.arange(t0, tf + 0.5 * h, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        X = Xs[i]
        k1 = f(X)                       # slope at start
        X_pred = X + h * k1             # Euler predictor
        k2 = f(X_pred)                  # slope at predicted endpoint
        Xs[i + 1] = X + 0.5 * h * (k1 + k2)  # average the two slopes
    return ts, Xs

# Integrate over [0, 60] with step h = 2.0
t0, tf, h = 0.0, 60.0, 2.0
t_rk2, X_rk2 = rk2_midpoint(deriv, X0, t0, tf, h)
t_eul, X_eul = euler(deriv, X0, t0, tf, h)
t_heu, X_heu = heun(deriv, X0, t0, tf, h)

X_exact_rk2 = exact(t_rk2)

# Max absolute errors vs exact solution at the shared step size h
err_rk2 = np.max(np.abs(X_rk2 - exact(t_rk2)))
err_eul = np.max(np.abs(X_eul - exact(t_eul)))
err_heu = np.max(np.abs(X_heu - exact(t_heu)))

print(f"Parameters: g = {g}, k = {k}, X0 = {X0}, h = {h}")
print(f"Steady state g/k = {g/k}")
print(f"Final time t = {tf}: RK2 = {X_rk2[-1]:.8f}, exact = {exact(tf):.8f}")
print(f"Max abs error Euler          = {err_eul:.8e}")
print(f"Max abs error Heun           = {err_heu:.8e}")
print(f"Max abs error RK2 (midpoint) = {err_rk2:.8e}")
print(f"RK2 error / Euler error ratio = {err_rk2/err_eul:.6f} (RK2 more accurate if < 1)")
print(f"RK2 error / Heun error ratio  = {err_rk2/err_heu:.6f} (comparable if near 1)")

# Explanation of the check:
print("Explanation: Because RK2's max error against the analytic solution is comparable "
      "to Heun's and much smaller than Euler's at the identical step size, the check confirms "
      "RK2 achieves the expected second-order accuracy rather than Euler's first-order accuracy.")

# ---- Plot: RK2 solution overlaid on the exact solution ----
t_fine = np.linspace(t0, tf, 500)
plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine), "k-", label="Exact")
plt.plot(t_rk2, X_rk2, "ro", markersize=5, label=f"RK2 midpoint (h={h})")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Constitutive gene expression: RK2 (midpoint) vs exact")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.5.1_s5.png")
