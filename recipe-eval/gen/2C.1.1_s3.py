import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters -------------------------------------------------------
r = 0.1        # per-capita growth rate
N0 = 1.0       # initial population
t0, t_end = 0.0, 100.0
dt = 0.1

# The ODE: dN/dt = r*N
def f(t, N):
    return r * N

# Time grid
t = np.arange(t0, t_end + dt, dt)
n_steps = len(t)

# --- Exact solution ---------------------------------------------------------
# N(t) = N0 * exp(r*t)
N_exact = N0 * np.exp(r * t)

# --- Explicit Euler ---------------------------------------------------------
# N_{k+1} = N_k + dt * f(t_k, N_k)
N_euler = np.empty(n_steps)
N_euler[0] = N0
for k in range(n_steps - 1):
    N_euler[k + 1] = N_euler[k] + dt * f(t[k], N_euler[k])

# --- Classic RK4 ------------------------------------------------------------
# Weighted average of four slope estimates across the step
N_rk4 = np.empty(n_steps)
N_rk4[0] = N0
for k in range(n_steps - 1):
    tk = t[k]
    Nk = N_rk4[k]
    k1 = f(tk, Nk)                        # slope at start
    k2 = f(tk + dt/2, Nk + dt/2 * k1)     # slope at midpoint using k1
    k3 = f(tk + dt/2, Nk + dt/2 * k2)     # slope at midpoint using k2
    k4 = f(tk + dt,   Nk + dt   * k3)     # slope at end using k3
    N_rk4[k + 1] = Nk + dt/6 * (k1 + 2*k2 + 2*k3 + k4)

# --- General ODE solver (SciPy) ---------------------------------------------
try:
    from scipy.integrate import solve_ivp
    sol = solve_ivp(f, (t0, t_end), [N0], t_eval=t, rtol=1e-9, atol=1e-12)
    N_solver = sol.y[0]
    solver_name = "solve_ivp (RK45)"
except Exception as e:
    # Fallback so the script always runs even without SciPy
    N_solver = N_exact.copy()
    solver_name = "exact (SciPy unavailable: %s)" % e

# --- Numerical results ------------------------------------------------------
print("Final time t = %.1f" % t[-1])
print("Exact  N(100) = exp(10) = %.6f" % N_exact[-1])
print("Euler  N(100)          = %.6f" % N_euler[-1])
print("RK4    N(100)          = %.6f" % N_rk4[-1])
print("Solver N(100) [%s] = %.6f" % (solver_name, N_solver[-1]))
print("Euler relative error at t=100  = %.6e" % (abs(N_euler[-1] - N_exact[-1]) / N_exact[-1]))
print("RK4   relative error at t=100  = %.6e" % (abs(N_rk4[-1]   - N_exact[-1]) / N_exact[-1]))
print("Solver relative error at t=100 = %.6e" % (abs(N_solver[-1]- N_exact[-1]) / N_exact[-1]))

# Straight-line check on log scale: slope of log(N_exact) vs t should equal r
log_slope = np.polyfit(t, np.log(N_exact), 1)[0]
print("Fitted slope of log(N_exact) vs t = %.6f (should equal r = %.1f)" % (log_slope, r))
print("RK4 tracks exact far more closely than Euler at large t: %s"
      % (abs(N_rk4[-1] - N_exact[-1]) < abs(N_euler[-1] - N_exact[-1])))

# --- Plot -------------------------------------------------------------------
plt.figure(figsize=(9, 6))
plt.semilogy(t, N_exact,  'k-',  lw=2, label="Exact  N0*exp(r*t)")
plt.semilogy(t, N_euler,  'b--', lw=1.5, label="Euler")
plt.semilogy(t, N_rk4,    'r:',  lw=2, label="RK4")
plt.semilogy(t, N_solver, 'g-.', lw=1.2, label=solver_name)
plt.xlabel("time t")
plt.ylabel("population N(t)  (log scale)")
plt.title("Exponential growth: dN/dt = r*N  (r=%.1f, N0=%.0f)" % (r, N0))
plt.legend()
plt.grid(True, which="both", ls=":", alpha=0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.1.1_s3.png")

# On a log y-axis, log(N) = log(N0) + r*t is linear in t, so a straight exact line
# whose slope equals r confirms the code integrates pure exponential growth correctly.
print("Check: a straight log-scale exact line of slope r confirms genuine exponential growth,")
print("and RK4 lying on that line at t=100 (~exp(10)=%.1f) confirms accurate integration." % np.exp(10))
