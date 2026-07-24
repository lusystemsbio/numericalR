import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Model: constitutively expressed gene
#   dX/dt = g - k*X   (transcription at rate g, degradation at rate k)
# ---------------------------------------------------------------
def f(t, X, g, k):
    return g - k * X

# Exact analytic solution for comparison / error measurement
def exact(t, g, k, X0):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ---------------------------------------------------------------
# Fourth-order Runge-Kutta, written out explicitly.
# Four slope evaluations per step: start, midpoint (twice), end,
# combined with weights 1, 2, 2, 1.
# ---------------------------------------------------------------
def rk4(f, t0, X0, h, n, g, k):
    t = np.empty(n + 1)
    X = np.empty(n + 1)
    t[0], X[0] = t0, X0
    for i in range(n):
        k1 = f(t[i],           X[i],                 g, k)  # slope at the start of the step
        k2 = f(t[i] + h/2,     X[i] + h/2 * k1,      g, k)  # slope at the midpoint using k1
        k3 = f(t[i] + h/2,     X[i] + h/2 * k2,      g, k)  # slope at the midpoint using k2
        k4 = f(t[i] + h,       X[i] + h   * k3,      g, k)  # slope at the end using k3
        # weighted average of the four slopes: (k1 + 2*k2 + 2*k3 + k4)/6
        X[i+1] = X[i] + h/6 * (k1 + 2*k2 + 2*k3 + k4)
        t[i+1] = t[i] + h
    return t, X

# ---------------------------------------------------------------
# Second-order Runge-Kutta (midpoint) for comparison
# ---------------------------------------------------------------
def rk2(f, t0, X0, h, n, g, k):
    t = np.empty(n + 1); X = np.empty(n + 1)
    t[0], X[0] = t0, X0
    for i in range(n):
        k1 = f(t[i],       X[i],            g, k)          # slope at the start
        k2 = f(t[i] + h/2, X[i] + h/2 * k1, g, k)          # slope at the midpoint
        X[i+1] = X[i] + h * k2                             # step using the midpoint slope
        t[i+1] = t[i] + h
    return t, X

# ---------------------------------------------------------------
# Forward Euler for comparison
# ---------------------------------------------------------------
def euler(f, t0, X0, h, n, g, k):
    t = np.empty(n + 1); X = np.empty(n + 1)
    t[0], X[0] = t0, X0
    for i in range(n):
        X[i+1] = X[i] + h * f(t[i], X[i], g, k)            # single slope at the start
        t[i+1] = t[i] + h
    return t, X

# ---------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------
g, k, X0 = 50.0, 0.1, 300.0
T = 60.0                      # total integration time
h = 2.0                       # step size for the main run
n = int(round(T / h))

# Integrate with RK4
t, X_rk4 = rk4(f, 0.0, X0, h, n, g, k)
X_true = exact(t, g, k, X0)

print(f"Steady state g/k = {g/k:.6f}")
print(f"RK4 final value X(T={T}) = {X_rk4[-1]:.6f}")
print(f"Exact final value X(T={T}) = {X_true[-1]:.6f}")

# ---------------------------------------------------------------
# Plot: RK4 solution overlaid on the exact solution
# ---------------------------------------------------------------
t_fine = np.linspace(0, T, 400)
plt.figure(figsize=(8, 5))
plt.plot(t_fine, exact(t_fine, g, k, X0), 'k-', lw=2, label="Exact solution")
plt.plot(t, X_rk4, 'ro', ms=6, label=f"RK4 (h={h})")
plt.xlabel("time t")
plt.ylabel("X(t)  (gene product)")
plt.title("RK4 vs exact solution:  dX/dt = g - k*X,  g=50, k=0.1, X0=300")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.6.1_s2.png")

# ---------------------------------------------------------------
# Accuracy check 1: RK4 much more accurate than RK2 and Euler
# at the SAME step size (max absolute error over the trajectory)
# ---------------------------------------------------------------
def max_err(method):
    _, Xn = method(f, 0.0, X0, h, n, g, k)
    return np.max(np.abs(Xn - exact(_, g, k, X0)))

err_euler = max_err(euler)
err_rk2   = max_err(rk2)
err_rk4   = max_err(rk4)
print(f"\nMax error at h={h}:")
print(f"  Euler max error = {err_euler:.6e}")
print(f"  RK2   max error = {err_rk2:.6e}")
print(f"  RK4   max error = {err_rk4:.6e}")

# ---------------------------------------------------------------
# Accuracy check 2: fourth-order convergence.
# Halving the step should cut the RK4 error by ~2^4 = 16.
# ---------------------------------------------------------------
def rk4_final_err(h_step):
    n_step = int(round(T / h_step))
    tt, XX = rk4(f, 0.0, X0, h_step, n_step, g, k)
    return abs(XX[-1] - exact(tt[-1], g, k, X0))

e_h  = rk4_final_err(h)
e_h2 = rk4_final_err(h / 2)
ratio = e_h / e_h2 if e_h2 != 0 else float('inf')
print(f"\nRK4 convergence (final-time error):")
print(f"  error at h={h}     = {e_h:.6e}")
print(f"  error at h/2={h/2} = {e_h2:.6e}")
print(f"  error ratio (h)/(h/2) = {ratio:.4f}  (expect ~16 for 4th order)")

# ---------------------------------------------------------------
# Explanation:
# The ratio near 16 = 2^4 confirms the error scales as h^4, which is the
# defining signature of a fourth-order method, so RK4 is behaving correctly.
# ---------------------------------------------------------------
print("\nExplanation: an error-reduction factor of ~16 when the step is halved")
print("means the global error scales as h^4, the hallmark of a genuine 4th-order")
print("method, which is exactly what a correct RK4 implementation must exhibit.")
