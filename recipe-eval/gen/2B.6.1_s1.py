import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model: constitutively expressed gene, dX/dt = g - k*X ----
def f(t, X, g, k):
    return g - k * X  # slope (RHS of the ODE)

# ---- Exact solution ----
def exact(t, g, k, X0):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ---- RK4 integrator (implemented explicitly) ----
def rk4(f, X0, t0, tf, h, g, k):
    n = int(round((tf - t0) / h))
    t = np.empty(n + 1)
    X = np.empty(n + 1)
    t[0], X[0] = t0, X0
    for i in range(n):
        k1 = f(t[i],           X[i],              g, k)  # slope at start
        k2 = f(t[i] + h/2,     X[i] + h/2 * k1,   g, k)  # slope at midpoint (using k1)
        k3 = f(t[i] + h/2,     X[i] + h/2 * k2,   g, k)  # slope at midpoint (using k2)
        k4 = f(t[i] + h,       X[i] + h   * k3,   g, k)  # slope at end
        # combine with weights 1, 2, 2, 1
        X[i+1] = X[i] + (h / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        t[i+1] = t[i] + h
    return t, X

# ---- RK2 (midpoint) integrator, for comparison ----
def rk2(f, X0, t0, tf, h, g, k):
    n = int(round((tf - t0) / h))
    t = np.empty(n + 1); X = np.empty(n + 1)
    t[0], X[0] = t0, X0
    for i in range(n):
        k1 = f(t[i],       X[i],            g, k)
        k2 = f(t[i] + h/2, X[i] + h/2 * k1, g, k)
        X[i+1] = X[i] + h * k2
        t[i+1] = t[i] + h
    return t, X

# ---- Euler integrator, for comparison ----
def euler(f, X0, t0, tf, h, g, k):
    n = int(round((tf - t0) / h))
    t = np.empty(n + 1); X = np.empty(n + 1)
    t[0], X[0] = t0, X0
    for i in range(n):
        X[i+1] = X[i] + h * f(t[i], X[i], g, k)
        t[i+1] = t[i] + h
    return t, X

# ---- Parameters ----
g, k, X0 = 50.0, 0.1, 300.0
t0, tf = 0.0, 60.0

# ---- Solve with RK4 and compare to exact ----
h = 2.0
t, Xrk4 = rk4(f, X0, t0, tf, h, g, k)
Xex = exact(t, g, k, X0)

# ---- Plot: RK4 overlaid on exact ----
tfine = np.linspace(t0, tf, 500)
plt.figure(figsize=(8, 5))
plt.plot(tfine, exact(tfine, g, k, X0), 'k-', label="Exact")
plt.plot(t, Xrk4, 'ro', markersize=6, label=f"RK4 (h={h})")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Constitutive gene expression: RK4 vs exact")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.6.1_s1.png")

# ---- Accuracy check: max error at same step size h ----
def max_err(method, h):
    tt, XX = method(f, X0, t0, tf, h, g, k)
    return np.max(np.abs(XX - exact(tt, g, k, X0)))

err_euler = max_err(euler, h)
err_rk2   = max_err(rk2,   h)
err_rk4   = max_err(rk4,   h)

print(f"Max error Euler (h={h}):     {err_euler:.6e}")
print(f"Max error RK2   (h={h}):     {err_rk2:.6e}")
print(f"Max error RK4   (h={h}):     {err_rk4:.6e}")
print(f"RK4 more accurate than RK2 by factor: {err_rk2/err_rk4:.3e}")
print(f"RK4 more accurate than Euler by factor: {err_euler/err_rk4:.3e}")

# ---- Order check: halving the step should cut RK4 error by ~16 (4th order) ----
err_rk4_h  = max_err(rk4, h)
err_rk4_h2 = max_err(rk4, h/2)
ratio = err_rk4_h / err_rk4_h2
print(f"Max error RK4 (h={h}):        {err_rk4_h:.6e}")
print(f"Max error RK4 (h={h/2}):      {err_rk4_h2:.6e}")
print(f"Error-reduction ratio when step halved: {ratio:.4f}")
print(f"Estimated order of accuracy (log2 of ratio): {np.log2(ratio):.4f}")

# Explanation: A ratio near 16 (= 2^4) confirms RK4 is fourth-order accurate,
# since a p-th order method's global error scales as h^p, so halving h reduces
# error by 2^p, and 2^4 = 16 is exactly what fourth-order predicts.
print("Check meaning: ratio ~16 = 2^4 confirms fourth-order accuracy, "
      "because halving h in an O(h^p) method reduces error by 2^p.")
