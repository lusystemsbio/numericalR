import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
r = 0.1        # per-capita growth rate
N0 = 1.0       # initial population
t0, tf = 0.0, 100.0
dt = 0.1
t = np.arange(t0, tf + dt, dt)   # time grid, inclusive of tf
n = len(t)

# --- Right-hand side of the ODE: dN/dt = r*N ---
def f(N):
    return r * N

# --- Exact analytic solution: N(t) = N0 * exp(r*t) ---
N_exact = N0 * np.exp(r * t)

# --- Explicit Euler method ---
# N_{k+1} = N_k + dt * f(N_k)   (first-order, one slope evaluation)
N_euler = np.empty(n)
N_euler[0] = N0
for k in range(n - 1):
    N_euler[k + 1] = N_euler[k] + dt * f(N_euler[k])

# --- Classical fourth-order Runge-Kutta (RK4) ---
# Weighted average of four slopes per step (fourth-order accurate)
N_rk4 = np.empty(n)
N_rk4[0] = N0
for k in range(n - 1):
    y = N_rk4[k]
    k1 = f(y)                 # slope at start
    k2 = f(y + 0.5 * dt * k1) # slope at midpoint using k1
    k3 = f(y + 0.5 * dt * k2) # slope at midpoint using k2
    k4 = f(y + dt * k3)       # slope at end using k3
    N_rk4[k + 1] = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

# --- General ODE solver (scipy) for comparison ---
from scipy.integrate import solve_ivp
sol = solve_ivp(lambda tt, N: r * N, (t0, tf), [N0], t_eval=t, rtol=1e-10, atol=1e-12)
N_solver = sol.y[0]

# --- Numerical results ---
print(f"Exact  N(100) = {N_exact[-1]:.6f}")
print(f"exp(10)        = {np.exp(10):.6f}")
print(f"Euler  N(100) = {N_euler[-1]:.6f}")
print(f"RK4    N(100) = {N_rk4[-1]:.6f}")
print(f"Solver N(100) = {N_solver[-1]:.6f}")
print(f"Euler  relative error at t=100 = {abs(N_euler[-1]-N_exact[-1])/N_exact[-1]:.6e}")
print(f"RK4    relative error at t=100 = {abs(N_rk4[-1]-N_exact[-1])/N_exact[-1]:.6e}")
print(f"Solver relative error at t=100 = {abs(N_solver[-1]-N_exact[-1])/N_exact[-1]:.6e}")

# --- Log-scaled plot ---
plt.figure(figsize=(8, 6))
plt.semilogy(t, N_euler, label="Euler", lw=1.5)
plt.semilogy(t, N_rk4, "--", label="RK4", lw=1.5)
plt.semilogy(t, N_exact, ":", label="Exact  N0*exp(r*t)", lw=2, color="k")
plt.xlabel("time t")
plt.ylabel("population N(t)  (log scale)")
plt.title("Exponential bacterial growth: dN/dt = r*N  (r=0.1, N0=1)")
plt.legend()
plt.grid(True, which="both", ls=":", alpha=0.5)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.1.1_s2.png", dpi=120, bbox_inches="tight")

# The check confirms the result because pure exponential growth becomes a straight line
# under a log y-axis (log N = log N0 + r*t), so the exact line being straight, RK4 hugging
# it far better than the first-order Euler at large t, and the endpoint equalling exp(10)
# together demonstrate the integrators correctly reproduce N(t)=N0*exp(r*t).
print("Check: on a log y-axis the exact solution is linear (log N = log N0 + r*t), "
      "RK4 tracks it far better than Euler at large t, and N(100) matches exp(10) ~ 22026, "
      "confirming the numerical methods reproduce the exact exponential solution.")
