import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model parameters -----
g = 50.0    # constant transcription rate (nM/time)
k = 0.1     # linear degradation rate (1/time)
X0 = 300.0  # initial expression level (nM)
t_end = 80.0

steady_state = g / k  # analytical steady state g/k

# ----- Exact solution: X(t) = g/k + (X0 - g/k)*exp(-k*t) -----
def exact(t):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ----- Explicit Euler integrator -----
def euler(dt):
    n = int(round(t_end / dt))          # number of steps
    ts = np.empty(n + 1)                # time points
    xs = np.empty(n + 1)                # numerical X values
    ts[0] = 0.0
    xs[0] = X0                          # start from initial value
    for i in range(n):
        # forward Euler step: X(t+dt) = X(t) + (g - k*X(t))*dt
        xs[i + 1] = xs[i] + (g - k * xs[i]) * dt
        ts[i + 1] = ts[i] + dt          # advance time
    return ts, xs

# Integrate at both step sizes
t1, x1 = euler(1.0)
t2, x2 = euler(0.1)

# Exact solution sampled on a fine grid (for plotting) and on Euler grids (for error)
t_fine = np.linspace(0, t_end, 1000)
x_fine = exact(t_fine)

err1 = np.max(np.abs(x1 - exact(t1)))   # max error at dt = 1
err2 = np.max(np.abs(x2 - exact(t2)))   # max error at dt = 0.1

# Gap near t = 10 for dt = 1 (visible discrepancy)
gap_at_10 = abs(x1[np.argmin(np.abs(t1 - 10.0))] - exact(10.0))

# ----- Report numerical results -----
print(f"Analytical steady state g/k = {steady_state:.4f} nM")
print(f"Euler final value (dt=1.0):   {x1[-1]:.6f} nM")
print(f"Euler final value (dt=0.1):   {x2[-1]:.6f} nM")
print(f"Exact final value (t=80):     {exact(t_end):.6f} nM")
print(f"Max absolute error (dt=1.0):  {err1:.6f} nM")
print(f"Max absolute error (dt=0.1):  {err2:.6f} nM")
print(f"Error at t=10 (dt=1.0):       {gap_at_10:.6f} nM")
print(f"Check dt=0.1 tracks exact (max err < 1 nM): {err2 < 1.0}")
print(f"Check dt=1.0 visible gap near t=10 (> 1 nM): {gap_at_10 > 1.0}")
print(f"Check both reach g/k=500 (within 5 nM): "
      f"{abs(x1[-1]-steady_state) < 5.0 and abs(x2[-1]-steady_state) < 5.0}")

# Explanation: A max error well under 1 nM at dt=0.1 alongside a visible gap at
# dt=1 with both curves converging to g/k=500 confirms the Euler code is correct
# because it reproduces the exact solution as the step size shrinks (the hallmark
# of a consistent, convergent method) while still showing the expected O(dt)
# truncation error at the coarse step.

# ----- Plot -----
plt.figure(figsize=(9, 6))
plt.plot(t_fine, x_fine, 'k-', lw=2, label='Exact')
plt.plot(t1, x1, 'ro--', ms=4, label='Euler dt = 1')
plt.plot(t2, x2, 'b-', lw=1, label='Euler dt = 0.1')
plt.axhline(steady_state, color='gray', ls=':', label='g/k = 500')
plt.xlabel('time')
plt.ylabel('X (nM)')
plt.title('Constitutive gene expression: Euler vs exact solution')
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.1.1_s3.png")
