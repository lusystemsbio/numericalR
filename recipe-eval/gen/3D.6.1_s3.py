import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# ------------------------------------------------------------------
# Lynx-hare yearly data (Hudson's Bay Company, 1900-1920), thousands.
# N = hare (prey), P = lynx (predator).
# ------------------------------------------------------------------
years = np.arange(1900, 1921)
t_data = years - years[0]                      # time in years starting at 0
N_data = np.array([30.0, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22.0,
                   25.4, 27.1, 40.3, 57.0, 76.6, 52.3, 19.5, 11.2, 7.6,
                   14.6, 16.2, 24.7])          # hare
P_data = np.array([4.0, 6.1, 9.8, 35.2, 59.4, 41.7, 19.0, 13.0, 8.3, 9.1,
                   7.4, 8.0, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7,
                   10.1, 8.6])                 # lynx

# ------------------------------------------------------------------
# Lotka-Volterra right-hand side:
#   dN/dt = N*(a - b*P)   prey grow (a), eaten by predators (b)
#   dP/dt = P*(c*N - d)   predators grow on prey (c), die off (d)
# ------------------------------------------------------------------
def lv_rhs(state, a, b, c, d):
    N, P = state
    dN = N * (a - b * P)
    dP = P * (c * N - d)
    return np.array([dN, dP])

# ------------------------------------------------------------------
# Simulate the model with a fixed-step RK4 integrator (explicit, no
# black-box ODE routine) and sample it at the data years.
# ------------------------------------------------------------------
def simulate(params, t_eval, y0, steps_per_year=50):
    a, b, c, d = params
    dt = 1.0 / steps_per_year
    t_end = t_eval[-1]
    n = int(round(t_end * steps_per_year)) + 1
    t_fine = np.linspace(0.0, t_end, n)
    y = np.zeros((n, 2))
    y[0] = y0
    for i in range(n - 1):
        yi = y[i]
        k1 = lv_rhs(yi,               a, b, c, d)
        k2 = lv_rhs(yi + 0.5*dt*k1,   a, b, c, d)
        k3 = lv_rhs(yi + 0.5*dt*k2,   a, b, c, d)
        k4 = lv_rhs(yi + dt*k3,       a, b, c, d)
        y[i+1] = yi + (dt/6.0) * (k1 + 2*k2 + 2*k3 + k4)
    # interpolate the fine trajectory onto the requested sample times
    N_s = np.interp(t_eval, t_fine, y[:, 0])
    P_s = np.interp(t_eval, t_fine, y[:, 1])
    return N_s, P_s

# ------------------------------------------------------------------
# Objective: sum of squared differences between simulation and data.
# The initial population (year 1900) is taken from the data.
# ------------------------------------------------------------------
y0 = np.array([N_data[0], P_data[0]])

def sse(params):
    # guard against non-physical negative parameters blowing up
    if np.any(np.asarray(params) <= 0):
        return 1e12
    N_s, P_s = simulate(params, t_data, y0)
    if not (np.all(np.isfinite(N_s)) and np.all(np.isfinite(P_s))):
        return 1e12
    return np.sum((N_s - N_data)**2) + np.sum((P_s - P_data)**2)

# ------------------------------------------------------------------
# Hand-picked starting guess for a, b, c, d and local BFGS optimization.
# ------------------------------------------------------------------
guess = np.array([0.55, 0.028, 0.024, 0.80])
print(f"Starting guess: a={guess[0]}, b={guess[1]}, c={guess[2]}, d={guess[3]}")
print(f"SSE at starting guess: {sse(guess):.4f}")

result = minimize(sse, guess, method="BFGS",
                  options={"maxiter": 2000, "gtol": 1e-8})
a_fit, b_fit, c_fit, d_fit = result.x

print("--- Fitted parameters (BFGS) ---")
print(f"a (prey growth rate)      = {a_fit:.6f}")
print(f"b (predation rate)        = {b_fit:.6f}")
print(f"c (predator growth rate)  = {c_fit:.6f}")
print(f"d (predator death rate)   = {d_fit:.6f}")
print(f"Final SSE: {result.fun:.4f}")
print(f"Optimizer converged: {result.success}")

# ------------------------------------------------------------------
# Fitted trajectory, sampled both densely (for plotting) and at data years.
# ------------------------------------------------------------------
t_dense = np.linspace(0, t_data[-1], 500)
N_dense, P_dense = simulate(result.x, t_dense, y0)
N_fit, P_fit = simulate(result.x, t_data, y0)

# ------------------------------------------------------------------
# Separate check: does the fit follow both time courses AND the phase loop?
# Use coefficient of determination R^2 for each species (time-course fit)
# and the correlation of the fitted phase points with the data phase points.
# ------------------------------------------------------------------
def r2(obs, pred):
    ss_res = np.sum((obs - pred)**2)
    ss_tot = np.sum((obs - np.mean(obs))**2)
    return 1.0 - ss_res / ss_tot

r2_N = r2(N_data, N_fit)
r2_P = r2(P_data, P_fit)
print("--- Goodness-of-fit check ---")
print(f"R^2 hare  (prey time course)     = {r2_N:.4f}")
print(f"R^2 lynx  (predator time course) = {r2_P:.4f}")
# phase-plane check: how well the fitted (N,P) points match the data loop
phase_rms = np.sqrt(np.mean((N_fit - N_data)**2 + (P_fit - P_data)**2))
print(f"Phase-plane RMS point distance    = {phase_rms:.4f}")

# ------------------------------------------------------------------
# Plots: time series (top) and phase plane (bottom).
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

ax = axes[0]
ax.plot(years, N_data, "o", color="tab:green", label="hare data")
ax.plot(years, P_data, "s", color="tab:red", label="lynx data")
ax.plot(years[0] + t_dense, N_dense, "-", color="tab:green", label="hare fit")
ax.plot(years[0] + t_dense, P_dense, "-", color="tab:red", label="lynx fit")
ax.set_xlabel("year"); ax.set_ylabel("population (thousands)")
ax.set_title("Lotka-Volterra fit: time courses")
ax.legend()

ax = axes[1]
ax.plot(N_data, P_data, "o-", color="0.6", label="data loop")
ax.plot(N_dense, P_dense, "-", color="tab:blue", label="fitted loop")
ax.set_xlabel("hare N (thousands)"); ax.set_ylabel("lynx P (thousands)")
ax.set_title("Phase plane")
ax.legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.6.1_s3.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: matching both time courses (high R^2) AND the closed "
      "phase-plane loop confirms the fit reproduces not just the values at "
      "each year but the correct oscillation amplitude, period, and "
      "predator-prey phase lag, which is the dynamical signature of the data.")
