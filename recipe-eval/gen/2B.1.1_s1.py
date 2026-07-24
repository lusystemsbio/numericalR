import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters (constitutively expressed gene)
g = 50.0    # constant transcription rate
k = 0.1     # linear degradation rate
X0 = 300.0  # initial expression level

t0, t_end = 0.0, 80.0
steady_state = g / k  # g/k = 500 nM

def euler(dt):
    # Build the time grid from t=0 to t_end with step dt
    n_steps = int(round((t_end - t0) / dt))
    t = t0 + dt * np.arange(n_steps + 1)
    X = np.empty(n_steps + 1)
    X[0] = X0  # initial value
    # Explicit Euler stepping: advance one dt at a time
    for i in range(n_steps):
        dXdt = g - k * X[i]              # right-hand side of the ODE at current state
        X[i + 1] = X[i] + dXdt * dt      # Euler update rule
    return t, X

def exact(t):
    # Analytic solution of dX/dt = g - k*X
    return g / k + (X0 - g / k) * np.exp(-k * t)

# Numerical solutions at the two step sizes
t1, X1 = euler(1.0)
t01, X01 = euler(0.1)

# Exact solutions evaluated on the same grids for error comparison
Xexact1 = exact(t1)
Xexact01 = exact(t01)

# Max absolute errors
max_err_dt1 = np.max(np.abs(X1 - Xexact1))
max_err_dt01 = np.max(np.abs(X01 - Xexact01))

# Gap near t = 10 for dt = 1
idx_t10_dt1 = np.argmin(np.abs(t1 - 10.0))
gap_t10_dt1 = abs(X1[idx_t10_dt1] - Xexact1[idx_t10_dt1])

# Final (near steady-state) values
final_dt1 = X1[-1]
final_dt01 = X01[-1]
final_exact = Xexact1[-1]

print(f"Steady state g/k: {steady_state:.6f} nM")
print(f"Max error dt=1.0: {max_err_dt1:.6f} nM")
print(f"Max error dt=0.1: {max_err_dt01:.6f} nM")
print(f"Gap near t=10 (dt=1.0): {gap_t10_dt1:.6f} nM")
print(f"Euler final value at t=80, dt=1.0: {final_dt1:.6f} nM")
print(f"Euler final value at t=80, dt=0.1: {final_dt01:.6f} nM")
print(f"Exact final value at t=80: {final_exact:.6f} nM")
print(f"dt=0.1 max error well under 1 nM: {max_err_dt01 < 1.0}")

# Plot
t_fine = np.linspace(t0, t_end, 1000)
plt.figure(figsize=(9, 6))
plt.plot(t_fine, exact(t_fine), 'k-', lw=2, label="Exact")
plt.plot(t1, X1, 'r--o', ms=4, label="Euler dt=1")
plt.plot(t01, X01, 'b-', lw=1, label="Euler dt=0.1")
plt.axhline(steady_state, color='gray', ls=':', label="g/k = 500")
plt.xlabel("Time t")
plt.ylabel("X (nM)")
plt.title("Gene-expression ODE: Euler vs exact solution")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.1.1_s1.png")

# The check confirms the result because a max error well under 1 nM at dt=0.1 shows the
# Euler steps converge to the analytic curve as dt shrinks, while the visible ~gap at
# dt=1 near the fast-changing region and both solutions reaching g/k=500 confirm the
# method is correct but only first-order accurate.
print("Check explanation: shrinking dt from 1 to 0.1 drives the max error well below 1 nM "
      "while both curves still settle at g/k=500, confirming the Euler scheme is correctly "
      "implemented and converges to the exact solution as the step size decreases.")
