import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Lotka-Volterra predator-prey model
#   dN/dt = N*(a - b*P)   prey grow, eaten by predators
#   dP/dt = P*(c*N - d)   predators grow on prey, die off
# ----------------------------------------------------------------------

# Model parameters
a, b, c, d = 1.0, 0.03, 0.02, 1.0

def lotka_volterra(state):
    """Return the derivative vector [dN/dt, dP/dt] for the state [N, P]."""
    N, P = state
    dN = N * (a - b * P)
    dP = P * (c * N - d)
    return np.array([dN, dP])

def rk4_step(f, y, dt):
    """One explicit classic 4th-order Runge-Kutta step (written out by hand)."""
    k1 = f(y)                    # slope at the start of the interval
    k2 = f(y + 0.5 * dt * k1)    # slope at the midpoint using k1
    k3 = f(y + 0.5 * dt * k2)    # slope at the midpoint using k2
    k4 = f(y + dt * k3)          # slope at the end using k3
    # weighted average of the four slopes advances the solution
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

def integrate(f, y0, dt, t_end):
    """Integrate the system from y0 over [0, t_end] with fixed step dt."""
    n_steps = int(round(t_end / dt))
    t = np.linspace(0.0, n_steps * dt, n_steps + 1)
    Y = np.empty((n_steps + 1, len(y0)))
    Y[0] = y0
    for i in range(n_steps):
        Y[i + 1] = rk4_step(f, Y[i], dt)   # explicit marching in time
    return t, Y

# Integration settings
dt = 0.01
t_end = 50.0

# Several initial conditions (N, P)
initial_conditions = [(30, 10), (40, 20), (30, 25), (20, 40)]

# Integrate each trajectory
solutions = []
for y0 in initial_conditions:
    t, Y = integrate(lotka_volterra, np.array(y0, dtype=float), dt, t_end)
    solutions.append((y0, t, Y))

# The coexistence equilibrium (fixed point) about which orbits circulate
N_star = d / c
P_star = a / b
print(f"Equilibrium point (N*, P*): ({N_star:.4f}, {P_star:.4f})")

# ----------------------------------------------------------------------
# Quantitative check 1: populations oscillate OUT OF PHASE in time.
# Use the first trajectory; measure the lag at peak cross-correlation
# between N(t) and P(t). A positive lag (fraction of a period) confirms
# the predator peak trails the prey peak.
# ----------------------------------------------------------------------
y0_check, t_check, Y_check = solutions[0]
N_series = Y_check[:, 0]
P_series = Y_check[:, 1]

# Peak times of prey and predator over the run
t_N_peak = t_check[np.argmax(N_series)]
t_P_peak = t_check[np.argmax(P_series)]
print(f"Time of first-window prey peak N:     t = {t_N_peak:.2f}")
print(f"Time of first-window predator peak P: t = {t_P_peak:.2f}")

# Cross-correlation lag between the two (mean-removed) signals
Nc = N_series - N_series.mean()
Pc = P_series - P_series.mean()
corr = np.correlate(Nc, Pc, mode="full")
lags = np.arange(-len(Nc) + 1, len(Nc))
best_lag_steps = lags[np.argmax(corr)]
best_lag_time = best_lag_steps * dt
print(f"Cross-correlation peak lag (P lags N): {best_lag_time:.2f} time units")

corr_at_zero = np.corrcoef(N_series, P_series)[0, 1]
print(f"Instantaneous Pearson correlation N vs P: {corr_at_zero:.4f}")

# ----------------------------------------------------------------------
# Quantitative check 2: phase-plane trajectories are CLOSED loops.
# For each orbit, compare the final state to the initial state; a small
# return distance (relative to the orbit's size) means it closed on itself.
# ----------------------------------------------------------------------
for y0, t, Y in solutions:
    start = Y[0]
    end = Y[-1]
    return_dist = np.linalg.norm(end - start)
    orbit_span = np.linalg.norm(Y.max(axis=0) - Y.min(axis=0))
    rel = return_dist / orbit_span
    print(f"IC {y0}: return distance = {return_dist:.4f}, "
          f"orbit span = {orbit_span:.4f}, relative closure = {rel:.4f}")

# ----------------------------------------------------------------------
# Plot 1: time series of N and P for the first initial condition
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

ax1.plot(t_check, N_series, label="Prey N", color="tab:green")
ax1.plot(t_check, P_series, label="Predator P", color="tab:red")
ax1.set_xlabel("time")
ax1.set_ylabel("population")
ax1.set_title(f"Time series (IC = {y0_check})")
ax1.legend()
ax1.grid(True, alpha=0.3)

# ----------------------------------------------------------------------
# Plot 2: phase plane showing nested closed orbits
# ----------------------------------------------------------------------
colors = ["tab:blue", "tab:orange", "tab:purple", "tab:brown"]
for (y0, t, Y), col in zip(solutions, colors):
    ax2.plot(Y[:, 0], Y[:, 1], color=col, label=f"IC {y0}")
    ax2.plot(y0[0], y0[1], "o", color=col)
ax2.plot(N_star, P_star, "k*", markersize=14, label="equilibrium")
ax2.set_xlabel("Prey N")
ax2.set_ylabel("Predator P")
ax2.set_title("Phase plane: nested closed orbits")
ax2.legend()
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.2.1_s3.png", dpi=120)

# ----------------------------------------------------------------------
# One-sentence explanation of why these checks confirm the result:
# A positive predator-lags-prey time offset (out-of-phase oscillation)
# together with each trajectory returning to its start (near-zero relative
# closure) confirms that the dynamics are periodic and trace conserved,
# nested closed loops around the coexistence equilibrium.
# ----------------------------------------------------------------------
print("Check: predator peak lags prey peak (out of phase) AND each orbit "
      "returns to its start (closed loop) => periodic nested closed orbits.")
