import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# -------------------------------------------------------------------
# Lynx-Hare yearly data (Hudson's Bay Co. pelt counts, thousands),
# 1900-1920. N = hare (prey), P = lynx (predator).
# -------------------------------------------------------------------
years = np.arange(1900, 1921)
t_data = years - years[0]                     # time in years starting at 0
N_data = np.array([30.0, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4,
                   22.0, 25.4, 27.1, 40.3, 57.0, 76.6, 52.3, 19.5,
                   11.2, 7.6, 14.6, 16.2, 24.7])   # hare
P_data = np.array([4.0, 6.1, 9.8, 35.2, 59.4, 41.7, 19.0, 13.0,
                   8.3, 9.1, 7.4, 8.0, 12.3, 19.5, 45.7, 51.1,
                   29.7, 15.8, 9.7, 10.1, 8.6])     # lynx

# -------------------------------------------------------------------
# Model: dN/dt = N*(a - b*P), dP/dt = P*(c*N - d)
# -------------------------------------------------------------------
def derivs(state, a, b, c, d):
    N, P = state
    dN = N * (a - b * P)
    dP = P * (c * N - d)
    return np.array([dN, dP])

# Explicit fixed-step RK4 integrator over the data times.
def simulate(params, t, y0, steps_per_year=20):
    a, b, c, d = params
    dt = (t[1] - t[0]) / steps_per_year        # sub-step size
    traj = np.zeros((len(t), 2))
    state = np.array(y0, dtype=float)
    traj[0] = state
    for i in range(1, len(t)):
        for _ in range(steps_per_year):        # step from t[i-1] to t[i]
            k1 = derivs(state, a, b, c, d)
            k2 = derivs(state + 0.5 * dt * k1, a, b, c, d)
            k3 = derivs(state + 0.5 * dt * k2, a, b, c, d)
            k4 = derivs(state + dt * k3, a, b, c, d)
            state = state + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            state = np.clip(state, 0.0, 1e6)   # keep populations non-negative/bounded
        traj[i] = state
    return traj

# -------------------------------------------------------------------
# Objective: sum of squared differences between simulation and data.
# -------------------------------------------------------------------
def sse(params):
    y0 = (N_data[0], P_data[0])                # start sim at first observation
    if np.any(np.array(params) <= 0):          # parameters must be positive
        return 1e12
    sim = simulate(params, t_data, y0)
    if not np.all(np.isfinite(sim)):
        return 1e12
    return np.sum((sim[:, 0] - N_data) ** 2) + np.sum((sim[:, 1] - P_data) ** 2)

# -------------------------------------------------------------------
# Hand-picked starting guess for (a, b, c, d) and BFGS optimization.
# -------------------------------------------------------------------
p0 = np.array([0.5, 0.02, 0.02, 0.8])         # reasonable initial guess
print(f"Starting guess: a={p0[0]}, b={p0[1]}, c={p0[2]}, d={p0[3]}")
print(f"Starting SSE: {sse(p0):.4f}")

result = minimize(sse, p0, method="BFGS",
                  options={"maxiter": 2000, "gtol": 1e-8})
a, b, c, d = result.x
print(f"Fitted a: {a:.6f}")
print(f"Fitted b: {b:.6f}")
print(f"Fitted c: {c:.6f}")
print(f"Fitted d: {d:.6f}")
print(f"Final SSE: {result.fun:.4f}")
print(f"Optimizer converged: {result.success}")

# -------------------------------------------------------------------
# Fitted trajectory on a fine time grid for smooth plotting.
# -------------------------------------------------------------------
t_fine = np.linspace(t_data[0], t_data[-1], 400)
sim_fine = simulate(result.x, t_fine, (N_data[0], P_data[0]),
                    steps_per_year=1)          # one RK4 step per fine interval
sim_at_data = simulate(result.x, t_data, (N_data[0], P_data[0]))

# -------------------------------------------------------------------
# Separate check: how well the fit follows the data, reported as the
# correlation between fitted and observed series (time course + loop).
# -------------------------------------------------------------------
def r2(obs, pred):
    ss_res = np.sum((obs - pred) ** 2)
    ss_tot = np.sum((obs - np.mean(obs)) ** 2)
    return 1.0 - ss_res / ss_tot

r2_hare = r2(N_data, sim_at_data[:, 0])
r2_lynx = r2(P_data, sim_at_data[:, 1])
print(f"R^2 hare (prey) time course: {r2_hare:.4f}")
print(f"R^2 lynx (predator) time course: {r2_lynx:.4f}")
print(f"Phase-plane check (corr hare vs lynx, data): "
      f"{np.corrcoef(N_data, P_data)[0,1]:.4f}")
print(f"Phase-plane check (corr hare vs lynx, fit):  "
      f"{np.corrcoef(sim_at_data[:,0], sim_at_data[:,1])[0,1]:.4f}")

# -------------------------------------------------------------------
# Plots: time courses and phase plane.
# -------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

ax = axes[0]
ax.plot(years, N_data, "o", color="tab:green", label="Hare data")
ax.plot(years, P_data, "s", color="tab:red", label="Lynx data")
ax.plot(years[0] + t_fine, sim_fine[:, 0], "-", color="tab:green",
        label="Hare fit")
ax.plot(years[0] + t_fine, sim_fine[:, 1], "-", color="tab:red",
        label="Lynx fit")
ax.set_xlabel("Year")
ax.set_ylabel("Population (thousands)")
ax.set_title("Lotka-Volterra fit: time courses")
ax.legend()

ax = axes[1]
ax.plot(N_data, P_data, "o-", color="0.5", alpha=0.7, label="Data loop")
ax.plot(sim_fine[:, 0], sim_fine[:, 1], "-", color="tab:blue",
        label="Fitted loop")
ax.set_xlabel("Hare (prey)")
ax.set_ylabel("Lynx (predator)")
ax.set_title("Phase plane")
ax.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.6.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Check rationale: if the fitted curves overlay the data in BOTH the "
      "time courses AND trace the same closed loop in the phase plane, then "
      "the model reproduces the timing, amplitude, and coupled predator-prey "
      "cycling of the observations, which is exactly what a correct fit must do.")
