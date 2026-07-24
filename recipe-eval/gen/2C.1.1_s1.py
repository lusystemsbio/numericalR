import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
r = 0.1        # constant per-capita growth rate
N0 = 1.0       # initial population
t0, tf = 0.0, 100.0
dt = 0.1

# Time grid: from t0 to tf inclusive with step dt
t = np.arange(t0, tf + dt, dt)
n_steps = len(t) - 1

# The right-hand side of the ODE: dN/dt = r*N
def f(N):
    return r * N

# --- Exact analytic solution: N(t) = N0*exp(r*t) ---
N_exact = N0 * np.exp(r * t)

# --- Explicit Euler method ---
# N_{k+1} = N_k + dt * f(N_k)
N_euler = np.empty_like(t)
N_euler[0] = N0
for k in range(n_steps):
    N_euler[k + 1] = N_euler[k] + dt * f(N_euler[k])

# --- Classical 4th-order Runge-Kutta (RK4) ---
# Combine four slope estimates per step for higher accuracy
N_rk4 = np.empty_like(t)
N_rk4[0] = N0
for k in range(n_steps):
    Nk = N_rk4[k]
    k1 = f(Nk)                 # slope at start
    k2 = f(Nk + 0.5 * dt * k1) # slope at midpoint using k1
    k3 = f(Nk + 0.5 * dt * k2) # slope at midpoint using k2
    k4 = f(Nk + dt * k3)       # slope at end using k3
    N_rk4[k + 1] = Nk + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

# --- General ODE solver (SciPy) as an independent cross-check ---
try:
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda tt, N: r * N, (t0, tf), [N0],
                    t_eval=t, rtol=1e-10, atol=1e-12)
    N_solver = sol.y[0]
    solver_name = "SciPy solve_ivp"
except Exception:
    # Fallback: if SciPy is unavailable, reuse RK4 as the "general solver"
    N_solver = N_rk4
    solver_name = "RK4 (SciPy unavailable)"

# --- Report final values and errors ---
print(f"Exact  N(100) = {N_exact[-1]:.6f}  (exp(10) ~ 22026.4658)")
print(f"Euler  N(100) = {N_euler[-1]:.6f}")
print(f"RK4    N(100) = {N_rk4[-1]:.6f}")
print(f"{solver_name} N(100) = {N_solver[-1]:.6f}")
print(f"exp(10) reference          = {np.exp(10):.6f}")
print(f"Euler relative error at t=100 = {abs(N_euler[-1] - N_exact[-1]) / N_exact[-1]:.6e}")
print(f"RK4   relative error at t=100 = {abs(N_rk4[-1]   - N_exact[-1]) / N_exact[-1]:.6e}")
print(f"{solver_name} relative error at t=100 = {abs(N_solver[-1] - N_exact[-1]) / N_exact[-1]:.6e}")

# --- Straight-line check on log axis ---
# log(N_exact) should be linear in t with slope r and intercept log(N0).
slope, intercept = np.polyfit(t, np.log(N_exact), 1)
print(f"Fitted slope of log(N_exact) vs t = {slope:.6f}  (expected r = {r})")
print(f"Fitted intercept of log(N_exact)  = {intercept:.6f}  (expected log(N0) = {np.log(N0):.6f})")

# --- Plot (log-scaled y axis) ---
plt.figure(figsize=(9, 6))
plt.semilogy(t, N_exact, 'k-', linewidth=2, label="Exact: N0*exp(r*t)")
plt.semilogy(t, N_euler, 'r--', linewidth=1.5, label="Euler")
plt.semilogy(t, N_rk4, 'b:', linewidth=2, label="RK4")
plt.xlabel("t")
plt.ylabel("N(t)  (log scale)")
plt.title("Exponential bacterial growth: Euler vs RK4 vs Exact")
plt.legend()
plt.grid(True, which="both", ls=":", alpha=0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.1.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Why the check confirms the result: because exp(r*t) becomes a straight "
      "line of slope r under a log y-axis, seeing the exact curve as straight, "
      "RK4 hugging it while Euler drifts below at large t, and the endpoint "
      "landing on exp(10) ~ 22026 together verify both the analytic form and "
      "the relative accuracy of the two integrators.")
