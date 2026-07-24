import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- Model: constitutively expressed gene ---
# dX/dt = g - k*X
g, k, X0 = 50.0, 0.1, 300.0

def f(X):
    return g - k * X  # RHS depends only on X here

def exact(t):
    # analytic solution of the linear ODE
    return g / k + (X0 - g / k) * np.exp(-k * t)

# --- Heun integrator (explicit, second order) ---
def heun(h, T):
    ts = np.arange(0.0, T + h / 2, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        s1 = f(Xs[i])                 # slope at start
        Xp = Xs[i] + h * s1           # Euler-predicted endpoint
        s2 = f(Xp)                    # slope at predicted endpoint
        Xs[i + 1] = Xs[i] + h * 0.5 * (s1 + s2)  # advance by average slope
    return ts, Xs

# --- Euler integrator (first order, for comparison) ---
def euler(h, T):
    ts = np.arange(0.0, T + h / 2, h)
    Xs = np.empty_like(ts)
    Xs[0] = X0
    for i in range(len(ts) - 1):
        Xs[i + 1] = Xs[i] + h * f(Xs[i])  # single starting slope
    return ts, Xs

T = 60.0

# --- Solve and plot ---
h = 2.0
t_heun, X_heun = heun(h, T)
t_exact = np.linspace(0.0, T, 500)
X_exact = exact(t_exact)

plt.figure(figsize=(8, 5))
plt.plot(t_exact, X_exact, 'k-', lw=2, label="exact")
plt.plot(t_heun, X_heun, 'ro--', ms=5, label=f"Heun (h={h})")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Constitutive gene expression: Heun vs exact")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.4.1_s4.png")

# --- Check 1: Heun beats Euler at the same step size ---
def max_abs_err(ts, Xs):
    return np.max(np.abs(Xs - exact(ts)))

_, Xe = euler(h, T)
err_euler = max_abs_err(t_heun, Xe)
err_heun = max_abs_err(t_heun, X_heun)

print(f"Max abs error, Euler (h={h}): {err_euler:.6e}")
print(f"Max abs error, Heun  (h={h}): {err_heun:.6e}")
print(f"Heun more accurate than Euler at same h: {err_heun < err_euler}")

# --- Check 2: halving h cuts Heun error by ~4x (second-order convergence) ---
th1, Xh1 = heun(h, T)
th2, Xh2 = heun(h / 2, T)
e1 = max_abs_err(th1, Xh1)
e2 = max_abs_err(th2, Xh2)
ratio = e1 / e2

print(f"Heun error at h   ={h}:   {e1:.6e}")
print(f"Heun error at h/2 ={h/2}: {e2:.6e}")
print(f"Error ratio (should be ~4): {ratio:.4f}")

# Explanation: a factor-of-4 error drop when the step is halved is the
# signature of O(h^2) global accuracy (2^2 = 4), confirming Heun is second order.
print("Explanation: halving h dropping the error ~4x means error scales as h^2, confirming Heun's second-order accuracy.")
