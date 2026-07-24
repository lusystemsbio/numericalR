import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
r = 0.1        # per-capita growth rate
N0 = 1.0       # initial population
t0, tf = 0.0, 100.0
dt = 0.1

# --- Time grid ---
t = np.arange(t0, tf + dt, dt)   # includes t = 100
n_steps = len(t)

# --- The right-hand side of the ODE: dN/dt = r*N ---
def f(N):
    return r * N

# --- Explicit Euler integration ---
# N_{k+1} = N_k + dt * f(N_k): step forward using the slope at the current point
N_euler = np.empty(n_steps)
N_euler[0] = N0
for k in range(n_steps - 1):
    N_euler[k + 1] = N_euler[k] + dt * f(N_euler[k])

# --- Classical 4th-order Runge-Kutta (RK4) integration ---
# Combine four slope samples (start, two midpoints, end) into a weighted average
N_rk4 = np.empty(n_steps)
N_rk4[0] = N0
for k in range(n_steps - 1):
    Nk = N_rk4[k]
    k1 = f(Nk)                    # slope at the start of the interval
    k2 = f(Nk + 0.5 * dt * k1)    # slope at the midpoint using k1
    k3 = f(Nk + 0.5 * dt * k2)    # slope at the midpoint using k2
    k4 = f(Nk + dt * k3)          # slope at the end using k3
    N_rk4[k + 1] = Nk + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

# --- General ODE solver (scipy) for comparison ---
try:
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda tt, N: r * N, [t0, tf], [N0], t_eval=t, rtol=1e-10, atol=1e-12)
    N_solver = sol.y[0]
    solver_name = "scipy solve_ivp (RK45)"
except Exception:
    # Fallback: reuse RK4 result if scipy is unavailable
    N_solver = N_rk4
    solver_name = "scipy unavailable; using RK4"

# --- Exact solution ---
N_exact = N0 * np.exp(r * t)

# --- Numerical results ---
print(f"Exact N(100) = exp(10) = {np.exp(10.0):.6f}")
print(f"Euler  final value N(100) = {N_euler[-1]:.6f}")
print(f"RK4    final value N(100) = {N_rk4[-1]:.6f}")
print(f"Solver final value N(100) [{solver_name}] = {N_solver[-1]:.6f}")
print(f"Euler  final relative error = {abs(N_euler[-1] - N_exact[-1]) / N_exact[-1]:.6e}")
print(f"RK4    final relative error = {abs(N_rk4[-1] - N_exact[-1]) / N_exact[-1]:.6e}")
print(f"Solver final relative error = {abs(N_solver[-1] - N_exact[-1]) / N_exact[-1]:.6e}")

# Check straight line on log axis: log(N_exact) should be linear in t with slope r
logN = np.log(N_exact)
slope = (logN[-1] - logN[0]) / (t[-1] - t[0])
print(f"Slope of log(N_exact) vs t = {slope:.6f} (should equal r = {r})")
print(f"Max deviation of log(N_exact) from straight line r*t + log(N0) = "
      f"{np.max(np.abs(logN - (r * t + np.log(N0)))):.6e}")

# --- Plot ---
plt.figure(figsize=(8, 6))
plt.semilogy(t, N_exact, 'k-', lw=2, label="Exact: N0*exp(r*t)")
plt.semilogy(t, N_euler, 'b--', lw=1.5, label="Euler")
plt.semilogy(t, N_rk4, 'r:', lw=2, label="RK4")
plt.xlabel("t")
plt.ylabel("N(t)  (log scale)")
plt.title("Exponential bacterial growth: dN/dt = r*N")
plt.legend()
plt.grid(True, which="both", ls=":", alpha=0.5)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.1.1_s4.png")

# Explanation: On a log y-axis the exact exponential appears as a straight line of slope r,
# so a method that stays on that line out to large t and lands on exp(10) at t=100 has
# faithfully reproduced the true exponential growth; RK4's tiny error versus Euler's visible
# droop below the line confirms RK4's far higher accuracy.
print("Check: a log-y plot turns N0*exp(r*t) into a straight line of slope r, so a scheme "
      "that stays on that line to t=100 and hits exp(10) has reproduced the exact solution; "
      "RK4 tracks it while Euler falls below, confirming RK4's superior accuracy.")
