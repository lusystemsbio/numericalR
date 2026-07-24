import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: constitutive gene expression, dX/dt = g - k*X ----
g = 50.0     # constant transcription rate
k = 0.1      # linear degradation rate
X0 = 300.0   # initial amount

def f(X):
    # right-hand side (slope) of the ODE
    return g - k * X

def exact(t):
    # analytic solution
    return g / k + (X0 - g / k) * np.exp(-k * t)

def euler(dt, T):
    # simple first-order explicit Euler for comparison
    n = int(round(T / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        X[i + 1] = X[i] + dt * f(X[i])   # advance using slope at start
    return t, X

def heun(dt, T):
    # second-order Heun (explicit trapezoidal / improved Euler)
    n = int(round(T / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    X = np.empty(n + 1)
    X[0] = X0
    for i in range(n):
        s1 = f(X[i])                     # slope at the current point
        X_pred = X[i] + dt * s1          # Euler step to a predicted endpoint
        s2 = f(X_pred)                    # slope evaluated at the predicted endpoint
        X[i + 1] = X[i] + dt * 0.5 * (s1 + s2)  # advance by the average of the two slopes
    return t, X

# ---- Run integrators ----
T = 50.0
dt = 2.0

t_h, X_h = heun(dt, T)
t_e, X_e = euler(dt, T)
X_true = exact(t_h)

# ---- Plot: Heun overlaid on exact ----
tt = np.linspace(0.0, T, 500)
plt.figure(figsize=(8, 5))
plt.plot(tt, exact(tt), 'k-', label='Exact')
plt.plot(t_h, X_h, 'ro--', label='Heun (dt=%.1f)' % dt)
plt.xlabel('t')
plt.ylabel('X(t)')
plt.title('Constitutive gene expression: Heun vs exact')
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.4.1_s2.png")

# ---- Check 1: Heun beats Euler at the same step size ----
err_heun = np.max(np.abs(X_h - X_true))
err_euler = np.max(np.abs(X_e - exact(t_e)))
print("Max error Euler (dt=%.1f):        %.6e" % (dt, err_euler))
print("Max error Heun  (dt=%.1f):        %.6e" % (dt, err_heun))
print("Heun more accurate than Euler:    %s" % (err_heun < err_euler))

# ---- Check 2: halving dt cuts Heun error ~4x (second-order convergence) ----
_, Xh_half = heun(dt / 2.0, T)
err_heun_half = np.max(np.abs(Xh_half - exact(np.linspace(0.0, T, len(Xh_half)))))
ratio = err_heun / err_heun_half
print("Max error Heun  (dt=%.1f):        %.6e" % (dt / 2.0, err_heun_half))
print("Error ratio (dt -> dt/2):         %.4f" % ratio)
print("Ratio near 4 confirms 2nd order:  %s" % (3.0 < ratio < 5.0))

# One sentence: an error-reduction factor of ~4 when the step is halved means the
# global error scales like dt^2, which is the defining signature of a second-order
# method, confirming Heun is genuinely order 2 (and thus more accurate than Euler).
print("Explanation: halving dt shrinking the error ~4x means error ~ dt^2, the hallmark of a second-order method.")
