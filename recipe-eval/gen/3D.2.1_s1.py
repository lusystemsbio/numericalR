import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Lotka-Volterra predator-prey model
#   dN/dt = N*(a - b*P)   prey grow, eaten by predators
#   dP/dt = P*(c*N - d)   predators grow on prey, die off
# ---------------------------------------------------------------

# Model parameters
a, b, c, d = 1.0, 0.03, 0.02, 1.0

def deriv(state, t):
    """Return the two-variable time derivative [dN/dt, dP/dt]."""
    N, P = state
    dN = N * (a - b * P)   # prey equation
    dP = P * (c * N - d)   # predator equation
    return np.array([dN, dP])

def rk4_step(f, y, t, dt):
    """One classic 4th-order Runge-Kutta step, written out explicitly."""
    k1 = f(y, t)                       # slope at start
    k2 = f(y + 0.5 * dt * k1, t + 0.5 * dt)  # slope at midpoint using k1
    k3 = f(y + 0.5 * dt * k2, t + 0.5 * dt)  # slope at midpoint using k2
    k4 = f(y + dt * k3, t + dt)        # slope at end using k3
    # weighted average of the four slopes
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)

def integrate(f, y0, dt, T):
    """Integrate the system from y0 over [0, T] using RK4."""
    n = int(round(T / dt))
    t = np.linspace(0.0, T, n + 1)
    ys = np.empty((n + 1, len(y0)))
    ys[0] = y0
    for i in range(n):
        ys[i + 1] = rk4_step(f, ys[i], t[i], dt)  # advance one step
    return t, ys

# Integration settings
dt = 0.01
T = 50.0
initial_conditions = [(30, 10), (40, 20), (30, 25), (20, 40)]

# Integrate every initial condition
solutions = []
for y0 in initial_conditions:
    t, ys = integrate(deriv, np.array(y0, dtype=float), dt, T)
    solutions.append((y0, t, ys))

# ---------------------------------------------------------------
# Plots
# ---------------------------------------------------------------
fig, (ax_ts, ax_pp) = plt.subplots(1, 2, figsize=(13, 5))

# Time-series plot (use the first initial condition)
y0, t, ys = solutions[0]
ax_ts.plot(t, ys[:, 0], label="Prey N")
ax_ts.plot(t, ys[:, 1], label="Predator P")
ax_ts.set_xlabel("time")
ax_ts.set_ylabel("population")
ax_ts.set_title(f"Time series, IC (N,P)={y0}")
ax_ts.legend()

# Phase-plane plot: all initial conditions -> nested closed loops
for y0, t, ys in solutions:
    ax_pp.plot(ys[:, 0], ys[:, 1], label=f"IC {y0}")
ax_pp.set_xlabel("Prey N")
ax_pp.set_ylabel("Predator P")
ax_pp.set_title("Phase plane: nested closed orbits")
ax_pp.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.2.1_s1.png")

# ---------------------------------------------------------------
# Checks
# ---------------------------------------------------------------
# Check 1: out-of-phase oscillation.
# Cross-correlate the (mean-removed) N and P series; the lag that maximizes
# correlation should be a nonzero fraction of the period -> a phase shift.
y0, t, ys = solutions[0]
N = ys[:, 0] - ys[:, 0].mean()
P = ys[:, 1] - ys[:, 1].mean()
corr = np.correlate(N, P, mode="full")
lags = np.arange(-len(N) + 1, len(N))
best_lag_steps = lags[np.argmax(corr)]
best_lag_time = best_lag_steps * dt

# Estimate prey period from peak-to-peak spacing of N
Nfull = ys[:, 0]
peaks = [i for i in range(1, len(Nfull) - 1)
         if Nfull[i] > Nfull[i - 1] and Nfull[i] > Nfull[i + 1]]
period = (t[peaks[1]] - t[peaks[0]]) if len(peaks) >= 2 else float("nan")
phase_fraction = best_lag_time / period if period == period else float("nan")

print(f"Prey oscillation period (approx): {period:.4f} time units")
print(f"Best N-vs-P correlation lag: {best_lag_time:.4f} time units")
print(f"Phase offset as fraction of period: {phase_fraction:.4f}")
print(f"Out-of-phase (nonzero lag)?: {abs(best_lag_time) > 0.1}")

# Check 2: closed orbits.
# For each trajectory, compare the final state to the initial state.
# A closed loop returns near its start after a whole number of periods.
print("Closed-orbit check (distance between start and nearest later return):")
for y0, t, ys in solutions:
    start = ys[0]
    # search after the first quarter to skip the immediate neighborhood
    tail = ys[len(ys) // 4:]
    dists = np.linalg.norm(tail - start, axis=1)
    min_return = dists.min()
    scale = np.linalg.norm(start)
    print(f"  IC {y0}: min return distance = {min_return:.4f} "
          f"({100 * min_return / scale:.2f}% of |start|)")

# Report the range of orbits to show they are nested (different amplitudes)
print("Orbit N-amplitude (max-min prey) per initial condition:")
for y0, t, ys in solutions:
    amp = ys[:, 0].max() - ys[:, 0].min()
    print(f"  IC {y0}: prey amplitude = {amp:.4f}")

# One-sentence explanation of why these checks confirm the result:
print("Explanation: a nonzero correlation lag between N and P confirms the "
      "populations peak at different times (out of phase), while each "
      "trajectory returning close to its starting point with a distinct "
      "amplitude confirms the phase-plane orbits are closed and nested.")
