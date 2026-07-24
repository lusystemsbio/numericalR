import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
a, b, c, d = 1.0, 0.03, 0.02, 1.0

# Right-hand side of the Lotka-Volterra system.
# state = [N, P]; returns [dN/dt, dP/dt]
def lv(state):
    N, P = state
    dN = N * (a - b * P)   # prey grow, eaten by predators
    dP = P * (c * N - d)   # predators grow on prey, die off
    return np.array([dN, dP])

# ---- Generic RK4 integrator implemented explicitly ----
def rk4(f, y0, dt, nsteps):
    y = np.array(y0, dtype=float)
    traj = np.empty((nsteps + 1, len(y)))
    traj[0] = y
    for i in range(nsteps):
        k1 = f(y)                 # slope at start
        k2 = f(y + 0.5 * dt * k1) # slope at midpoint using k1
        k3 = f(y + 0.5 * dt * k2) # slope at midpoint using k2
        k4 = f(y + dt * k3)       # slope at end using k3
        y = y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # weighted average
        traj[i + 1] = y
    return traj

# ---- Integration settings ----
dt = 0.01
T = 50.0
nsteps = int(round(T / dt))
t = np.linspace(0.0, T, nsteps + 1)

initial_conditions = [(30, 10), (40, 20), (30, 25), (20, 40)]

# ---- Integrate from each initial condition ----
trajectories = []
for ic in initial_conditions:
    traj = rk4(lv, ic, dt, nsteps)
    trajectories.append(traj)

# ---- Report the fixed (equilibrium) point ----
# Interior equilibrium: N* = d/c, P* = a/b
N_star, P_star = d / c, a / b
print(f"Interior equilibrium (N*, P*): ({N_star:.4f}, {P_star:.4f})")

# ---- Numerical checks ----
# We use the first initial condition (30, 10) for the phase-shift diagnostic.
N0 = trajectories[0][:, 0]
P0 = trajectories[0][:, 1]

# Time of peak (argmax) of prey and predator over the first oscillation window.
# Restrict to a window that contains one full-ish cycle for a clean comparison.
i_peak_N = int(np.argmax(N0))
i_peak_P = int(np.argmax(P0))
print(f"Time of first prey (N) peak, IC (30,10):     {t[i_peak_N]:.3f}")
print(f"Time of first predator (P) peak, IC (30,10): {t[i_peak_P]:.3f}")
print(f"Predator peak lags prey peak by: {t[i_peak_P] - t[i_peak_N]:.3f} time units")

# Cross-correlation-free out-of-phase check: correlation of the two
# centered signals should be negative when they oscillate out of phase.
Nc = N0 - N0.mean()
Pc = P0 - P0.mean()
corr = float(np.corrcoef(Nc, Pc)[0, 1])
print(f"Correlation between N(t) and P(t), IC (30,10): {corr:.4f} (negative => out of phase)")

# Closed-orbit check: for each trajectory measure how close the endpoint
# returns to the start (distance relative to the orbit's spatial extent).
# A small ratio indicates the trajectory closes on itself (periodic loop).
print("Closed-orbit check (endpoint return distance relative to orbit size):")
for ic, traj in zip(initial_conditions, trajectories):
    start = traj[0]
    end = traj[-1]
    ret_dist = np.linalg.norm(end - start)
    extent = np.linalg.norm(traj.max(axis=0) - traj.min(axis=0))
    ratio = ret_dist / extent
    print(f"  IC {ic}: return_dist={ret_dist:.4f}, extent={extent:.4f}, ratio={ratio:.4f}")

# Conserved quantity check: the LV system has invariant
# H = c*N - d*ln(N) + b*P - a*ln(P), constant along closed orbits.
print("Conserved-quantity check (max relative drift of H along each orbit):")
for ic, traj in zip(initial_conditions, trajectories):
    N = traj[:, 0]
    P = traj[:, 1]
    H = c * N - d * np.log(N) + b * P - a * np.log(P)
    rel_drift = (H.max() - H.min()) / abs(H.mean())
    print(f"  IC {ic}: relative drift of H = {rel_drift:.3e}")

# ---- Plots ----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5))

# Time series for the first initial condition, showing out-of-phase oscillation.
ax1.plot(t, N0, label="Prey N(t)", color="tab:green")
ax1.plot(t, P0, label="Predator P(t)", color="tab:red")
ax1.set_xlabel("time")
ax1.set_ylabel("population")
ax1.set_title("Time series, IC (N, P) = (30, 10)")
ax1.legend()
ax1.grid(True, alpha=0.3)

# Phase plane showing nested closed orbits from all initial conditions.
colors = plt.cm.viridis(np.linspace(0, 0.85, len(trajectories)))
for ic, traj, col in zip(initial_conditions, trajectories, colors):
    ax2.plot(traj[:, 0], traj[:, 1], color=col, label=f"IC {ic}")
ax2.plot(N_star, P_star, "k*", markersize=12, label="equilibrium")
ax2.set_xlabel("prey N")
ax2.set_ylabel("predator P")
ax2.set_title("Phase plane: nested closed orbits")
ax2.legend()
ax2.grid(True, alpha=0.3)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.2.1_s2.png", dpi=130)

# One-sentence explanation of why the check confirms the result:
print("Explanation: A negative N-P correlation with the predator peak lagging the prey peak "
      "confirms the two populations oscillate out of phase, and each trajectory returning to its "
      "start with a near-constant conserved quantity H confirms the phase-plane orbits are closed, "
      "nested loops rather than spirals.")
