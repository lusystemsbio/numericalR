import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters: constitutive gene expression dX/dt = g - k*X
g = 50.0      # transcription rate (nM/time)
k = 0.1       # linear degradation rate (1/time)
X0 = 300.0    # initial expression level (nM)
t_start, t_end = 0.0, 80.0

def exact(t):
    # Analytic solution X(t) = g/k + (X0 - g/k)*exp(-k*t)
    return g / k + (X0 - g / k) * np.exp(-k * t)

def euler(dt):
    # Explicit Euler: step from t=0 to t_end with fixed dt
    n = int(round((t_end - t_start) / dt))   # number of steps
    t = np.empty(n + 1)
    X = np.empty(n + 1)
    t[0], X[0] = t_start, X0
    for i in range(n):
        # dX/dt evaluated at the current state
        deriv = g - k * X[i]
        # advance one step: X(t+dt) = X(t) + deriv*dt
        X[i + 1] = X[i] + deriv * dt
        t[i + 1] = t[i] + dt
    return t, X

# Run Euler at both step sizes
t1, X1 = euler(1.0)     # coarse step
t01, X01 = euler(0.1)   # fine step

# Steady state and max errors vs exact solution
steady = g / k
err1 = np.max(np.abs(X1 - exact(t1)))
err01 = np.max(np.abs(X01 - exact(t01)))

print(f"Steady state g/k = {steady:.4f} nM")
print(f"Euler dt=1   final X(80) = {X1[-1]:.6f} nM")
print(f"Euler dt=0.1 final X(80) = {X01[-1]:.6f} nM")
print(f"Exact        X(80)       = {exact(t_end):.6f} nM")
print(f"Max abs error, dt=1   : {err1:.6f} nM")
print(f"Max abs error, dt=0.1 : {err01:.6f} nM")
print(f"dt=0.1 max error well under 1 nM: {err01 < 1.0}")

# Plot: Euler solutions overlaid on exact
t_fine = np.linspace(t_start, t_end, 1000)
plt.figure(figsize=(9, 5.5))
plt.plot(t_fine, exact(t_fine), 'k-', lw=2.5, label="Exact")
plt.plot(t1, X1, 'ro-', ms=4, lw=1, label="Euler dt=1")
plt.plot(t01, X01, 'b--', lw=1.5, label="Euler dt=0.1")
plt.axhline(steady, color='gray', ls=':', label="g/k = 500")
plt.xlabel("time")
plt.ylabel("X (nM)")
plt.title("Constitutive gene expression: Euler vs exact")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.1.1_s2.png")

# One-sentence explanation of the check:
print("Check rationale: a max error well under 1 nM at dt=0.1 (vs a visible "
      "gap at dt=1) with both curves converging to g/k=500 confirms the Euler "
      "code is correct and that its error shrinks as the step size decreases, "
      "as expected for a first-order method.")
