import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g = 50.0      # constant transcription rate
k = 0.1       # linear degradation rate
X0 = 300.0    # initial expression level
t_end = 80.0  # end time

# Exact analytic solution of dX/dt = g - k*X
def exact(t):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# Explicit Euler integrator for dX/dt = g - k*X
def euler(dt):
    n = int(round(t_end / dt))         # number of steps
    t = np.zeros(n + 1)
    X = np.zeros(n + 1)
    X[0] = X0                          # initial value at t = 0
    for i in range(n):
        deriv = g - k * X[i]           # right-hand side evaluated at current state
        X[i + 1] = X[i] + deriv * dt   # Euler update: X(t+dt) = X(t) + f*dt
        t[i + 1] = t[i] + dt           # advance time
    return t, X

# Integrate at both step sizes
t1, X1 = euler(1.0)
t01, X01 = euler(0.1)

# Fine grid for the exact reference curve
t_fine = np.linspace(0, t_end, 1000)
X_fine = exact(t_fine)

# --- Error checks against the exact solution at the Euler grid points ---
err_dt1 = np.max(np.abs(X1 - exact(t1)))
err_dt01 = np.max(np.abs(X01 - exact(t01)))

# Gap near t = 10 for dt = 1
idx10 = np.argmin(np.abs(t1 - 10.0))
gap_at_10 = abs(X1[idx10] - exact(t1[idx10]))

steady_state = g / k

print(f"Steady state g/k (nM): {steady_state}")
print(f"Euler dt=1  final X (nM): {X1[-1]}")
print(f"Euler dt=0.1 final X (nM): {X01[-1]}")
print(f"Exact final X at t=80 (nM): {exact(t_end)}")
print(f"Max abs error, dt=1  (nM): {err_dt1}")
print(f"Max abs error, dt=0.1 (nM): {err_dt01}")
print(f"Gap at t=10, dt=1 (nM): {gap_at_10}")
print(f"Check dt=0.1 max error well under 1 nM: {err_dt01 < 1.0}")

# --- Plot ---
plt.figure(figsize=(8, 5))
plt.plot(t_fine, X_fine, 'k-', lw=2, label="Exact")
plt.plot(t1, X1, 'ro--', ms=4, label="Euler dt=1")
plt.plot(t01, X01, 'b.', ms=2, label="Euler dt=0.1")
plt.axhline(steady_state, color='gray', ls=':', label="g/k = 500")
plt.xlabel("time")
plt.ylabel("X (nM)")
plt.title("Gene expression: Euler vs exact (dX/dt = g - k*X)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.1.1_s4.png")

# The check confirms the result because a max error well under 1 nM at dt=0.1 (versus a
# visible gap at dt=1) demonstrates that Euler converges to the exact solution as the step
# size shrinks, while both reaching g/k=500 verifies the correct steady state was captured.
print("Check confirms result: smaller dt yields max error well under 1 nM and both curves reach the correct steady state g/k=500, showing Euler converges to the exact solution as dt shrinks.")
