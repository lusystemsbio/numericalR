import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
r = 0.1          # per-capita growth rate
N0 = 1.0         # initial population
t0, tf = 0.0, 100.0
dt = 0.1         # time step

# Time grid
t = np.arange(t0, tf + dt, dt)
n_steps = len(t)

# ---- Exact solution: N(t) = N0 * exp(r*t) ----
N_exact = N0 * np.exp(r * t)

# ---- Euler method (explicit, one step at a time) ----
# Update rule: N_{k+1} = N_k + dt * f(N_k), where f(N) = r*N
N_euler = np.empty(n_steps)
N_euler[0] = N0
for k in range(n_steps - 1):
    dNdt = r * N_euler[k]              # slope at current point
    N_euler[k + 1] = N_euler[k] + dt * dNdt

# ---- RK4 method (explicit, four-stage Runge-Kutta) ----
def f(N):
    return r * N                       # right-hand side dN/dt = r*N

N_rk4 = np.empty(n_steps)
N_rk4[0] = N0
for k in range(n_steps - 1):
    Nk = N_rk4[k]
    k1 = f(Nk)                         # slope at start
    k2 = f(Nk + 0.5 * dt * k1)         # slope at midpoint using k1
    k3 = f(Nk + 0.5 * dt * k2)         # slope at midpoint using k2
    k4 = f(Nk + dt * k3)               # slope at end using k3
    N_rk4[k + 1] = Nk + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # weighted avg

# ---- General ODE solver (SciPy solve_ivp, RK45) ----
from scipy.integrate import solve_ivp
sol = solve_ivp(lambda tt, N: r * N, (t0, tf), [N0], t_eval=t, rtol=1e-10, atol=1e-12)
N_solver = sol.y[0]

# ---- Log-scaled plot: Euler, RK4, exact ----
plt.figure(figsize=(9, 6))
plt.semilogy(t, N_exact, 'k-', lw=2, label="Exact: N0*exp(r*t)")
plt.semilogy(t, N_euler, 'r--', lw=1.5, label="Euler")
plt.semilogy(t, N_rk4, 'b:', lw=1.8, label="RK4")
plt.semilogy(t, N_solver, 'g-.', lw=1.2, label="solve_ivp (RK45)")
plt.xlabel("time t")
plt.ylabel("population N(t)  (log scale)")
plt.title("Exponential bacterial growth: dN/dt = r*N")
plt.legend()
plt.grid(True, which="both", ls=":", alpha=0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.1.1_s5.png")

# ---- Checks / numerical results ----
# Straightness of exact solution on log axis: log(N) = log(N0) + r*t is linear.
# Fit log(N_exact) vs t and report slope (should equal r) and R^2 (should be 1).
logN = np.log(N_exact)
slope, intercept = np.polyfit(t, logN, 1)
resid = logN - (slope * t + intercept)
ss_res = np.sum(resid**2)
ss_tot = np.sum((logN - logN.mean())**2)
r2 = 1.0 - ss_res / ss_tot

N100_theory = np.exp(10.0)

print(f"Fitted slope of log(N_exact) vs t (should equal r={r}): {slope:.12f}")
print(f"R^2 of linear fit on log axis (1.0 => straight line): {r2:.15f}")
print(f"Theoretical N(100) = exp(10): {N100_theory:.6f}")
print(f"Exact    N(100): {N_exact[-1]:.6f}")
print(f"Euler    N(100): {N_euler[-1]:.6f}")
print(f"RK4      N(100): {N_rk4[-1]:.6f}")
print(f"solve_ivp N(100): {N_solver[-1]:.6f}")
print(f"Euler relative error at t=100: {abs(N_euler[-1]-N_exact[-1])/N_exact[-1]:.6e}")
print(f"RK4   relative error at t=100: {abs(N_rk4[-1]-N_exact[-1])/N_exact[-1]:.6e}")
print(f"solve_ivp relative error at t=100: {abs(N_solver[-1]-N_exact[-1])/N_exact[-1]:.6e}")

# Explanation:
# The check confirms the result because on a log y-axis the exact exponential
# becomes the straight line log(N)=log(N0)+r*t (slope r, R^2=1), and RK4's
# far smaller error than Euler plus agreement with N(100)=exp(10)~22026 shows
# the explicit integrators are converging to the true analytic solution.
print("Why: on a log axis the exact solution is exactly linear (slope r, R^2=1), "
      "and RK4 hugging that line and hitting exp(10) while Euler drifts confirms "
      "the numerics reproduce the analytic exponential growth.")
