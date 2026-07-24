import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
g = 50.0     # constant transcription rate
k = 0.1      # linear degradation rate
X0 = 300.0   # initial value
t_end = 80.0 # end time

steady_state = g / k  # analytic steady state g/k
print(f"Steady state g/k: {steady_state:.6f} nM")

# ---- Exact solution: X(t) = g/k + (X0 - g/k)*exp(-k*t) ----
def exact(t):
    return g / k + (X0 - g / k) * np.exp(-k * t)

# ---- Explicit Euler integrator ----
def euler(dt):
    n = int(round(t_end / dt))          # number of steps
    t = np.empty(n + 1)                 # time array
    X = np.empty(n + 1)                 # solution array
    t[0], X[0] = 0.0, X0                # initial condition at t = 0
    for i in range(n):
        dXdt = g - k * X[i]             # right-hand side dX/dt = g - k*X
        X[i + 1] = X[i] + dXdt * dt     # Euler step: X(t+dt) = X(t) + f*dt
        t[i + 1] = t[i] + dt            # advance time
    return t, X

# ---- Integrate at both step sizes ----
t1, X1 = euler(1.0)     # dt = 1
t01, X01 = euler(0.1)   # dt = 0.1

# ---- Fine grid for the exact curve ----
t_fine = np.linspace(0, t_end, 2000)
X_exact_fine = exact(t_fine)

# ---- Error checks: compare Euler to exact at the same time points ----
err_dt1 = np.max(np.abs(X1 - exact(t1)))
err_dt01 = np.max(np.abs(X01 - exact(t01)))
print(f"Max error, dt = 1.0:  {err_dt1:.6f} nM")
print(f"Max error, dt = 0.1:  {err_dt01:.6f} nM")

# Gap near t = 10 for dt = 1
i10 = np.argmin(np.abs(t1 - 10.0))
gap_at_t10 = abs(X1[i10] - exact(t1[i10]))
print(f"Gap near t = 10 (dt = 1): {gap_at_t10:.6f} nM")

# Final values (both should approach g/k = 500)
print(f"Final value, dt = 1.0:  {X1[-1]:.6f} nM")
print(f"Final value, dt = 0.1:  {X01[-1]:.6f} nM")
print(f"Final value, exact:     {exact(t_end):.6f} nM")

print(f"dt = 0.1 error well under 1 nM: {err_dt01 < 1.0}")

# ---- Plot ----
plt.figure(figsize=(9, 6))
plt.plot(t_fine, X_exact_fine, 'k-', lw=2, label='Exact')
plt.plot(t1, X1, 'o-', ms=4, color='tab:red', label='Euler, dt = 1')
plt.plot(t01, X01, '-', color='tab:blue', label='Euler, dt = 0.1')
plt.axhline(steady_state, color='gray', ls='--', lw=1, label='g/k = 500')
plt.xlabel('time t')
plt.ylabel('X (nM)')
plt.title('Constitutive gene expression: Euler vs exact solution')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.1.1_s5.png")

# The check confirms the result because a max error well under 1 nM at dt = 0.1
# (versus a visible gap at dt = 1) shows Euler's error shrinks with step size and
# converges to the exact solution, while both trajectories relaxing to g/k = 500
# confirms the correct steady state is reached.
print("Check confirms result: smaller dt gives max error < 1 nM (vs visible gap at dt=1), and both reach g/k=500, demonstrating convergence to the exact solution.")
