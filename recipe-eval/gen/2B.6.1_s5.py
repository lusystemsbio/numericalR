import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: constitutively expressed gene ---
# dX/dt = g - k*X  (transcription at constant rate g, linear degradation at rate k)
g = 50.0
k = 0.1
X0 = 300.0

def f(t, X):
    # RHS of the ODE (does not depend on t here, but kept general)
    return g - k * X

def exact(t):
    # Analytic solution X(t) = g/k + (X0 - g/k)*exp(-k*t)
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Explicit fourth-order Runge-Kutta integrator ---
def rk4(f, X0, t0, tf, h):
    n = int(round((tf - t0) / h))
    ts = np.empty(n + 1)
    Xs = np.empty(n + 1)
    ts[0], Xs[0] = t0, X0
    t, X = t0, X0
    for i in range(n):
        k1 = f(t, X)                      # slope at the start of the step
        k2 = f(t + h / 2, X + h / 2 * k1) # slope at the midpoint using k1
        k3 = f(t + h / 2, X + h / 2 * k2) # slope at the midpoint using k2
        k4 = f(t + h, X + h * k3)         # slope at the end of the step
        # combine with weights 1, 2, 2, 1 (divided by 6)
        X = X + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        t = t + h
        ts[i + 1], Xs[i + 1] = t, X
    return ts, Xs

# --- Second-order Runge-Kutta (midpoint) for comparison ---
def rk2(f, X0, t0, tf, h):
    n = int(round((tf - t0) / h))
    t, X = t0, X0
    for i in range(n):
        k1 = f(t, X)                       # slope at start
        k2 = f(t + h / 2, X + h / 2 * k1)  # slope at midpoint
        X = X + h * k2                     # step using midpoint slope
        t = t + h
    return X

# --- Forward Euler for comparison ---
def euler(f, X0, t0, tf, h):
    n = int(round((tf - t0) / h))
    t, X = t0, X0
    for i in range(n):
        X = X + h * f(t, X)  # single slope at start
        t = t + h
    return X

# --- Integrate and plot ---
t0, tf, h = 0.0, 60.0, 2.0
ts, Xs = rk4(f, X0, t0, tf, h)
Xe = exact(ts)

plt.figure(figsize=(8, 5))
plt.plot(ts, Xe, "b-", lw=2, label="Exact")
plt.plot(ts, Xs, "ro", ms=5, label="RK4 (h=%.1f)" % h)
plt.xlabel("time")
plt.ylabel("X(t)  (gene product)")
plt.title("Constitutive gene expression: RK4 vs exact")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.6.1_s5.png")

# --- Accuracy checks (final-time error at t = tf) ---
Xt = exact(tf)
err_euler = abs(euler(f, X0, t0, tf, h) - Xt)
err_rk2 = abs(rk2(f, X0, t0, tf, h) - Xt)
err_rk4 = abs(rk4(f, X0, t0, tf, h)[1][-1] - Xt)

# RK4 error at h and at h/2, to check ~16x error reduction (4th order => 2^4)
err_rk4_h = abs(rk4(f, X0, t0, tf, h)[1][-1] - Xt)
err_rk4_h2 = abs(rk4(f, X0, t0, tf, h / 2)[1][-1] - Xt)
ratio = err_rk4_h / err_rk4_h2 if err_rk4_h2 != 0 else float("inf")

print("Exact final value X(%.1f)      = %.10f" % (tf, Xt))
print("RK4 final value                = %.10f" % Xs[-1])
print("Euler error at h=%.1f          = %.6e" % (h, err_euler))
print("RK2 error at h=%.1f            = %.6e" % (h, err_rk2))
print("RK4 error at h=%.1f            = %.6e" % (h, err_rk4))
print("RK4 error at h/2=%.1f          = %.6e" % (h / 2, err_rk4_h2))
print("Error reduction ratio (h -> h/2) = %.4f  (expect ~16 for 4th order)" % ratio)

# One-sentence explanation:
# Because RK4 is a fourth-order method, its global error scales like h^4, so
# halving the step multiplies the error by (1/2)^4 = 1/16; observing the error
# drop by ~16x (and being far smaller than RK2's ~h^2 and Euler's ~h^1 errors)
# confirms the integrator is behaving as a correct fourth-order scheme.
print("Explanation: RK4 error ~ h^4, so halving h cuts error by 2^4=16, and its")
print("far smaller error than RK2 (~h^2) and Euler (~h^1) confirms 4th-order accuracy.")
