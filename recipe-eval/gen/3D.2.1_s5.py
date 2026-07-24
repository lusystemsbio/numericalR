import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
a, b, c, d = 1.0, 0.03, 0.02, 1.0

# Right-hand side of the Lotka-Volterra system.
# state = [N (prey), P (predators)]
def lv(state):
    N, P = state
    dN = N * (a - b * P)   # prey grow, eaten by predators
    dP = P * (c * N - d)   # predators grow on prey, die off
    return np.array([dN, dP])

# --- Generic RK4 integrator (implemented explicitly) ---
def rk4_step(f, y, dt):
    k1 = f(y)                 # slope at start
    k2 = f(y + 0.5 * dt * k1) # slope at midpoint using k1
    k3 = f(y + 0.5 * dt * k2) # slope at midpoint using k2
    k4 = f(y + dt * k3)       # slope at end using k3
    return y + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # weighted average

def integrate(f, y0, dt, T):
    n = int(round(T / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    Y = np.empty((n + 1, len(y0)))
    Y[0] = y0
    for i in range(n):
        Y[i + 1] = rk4_step(f, Y[i], dt)
    return t, Y

# --- Simulation settings ---
dt, T = 0.01, 50.0
initial_conditions = [(30, 10), (40, 20), (30, 25), (20, 40)]

results = []
for ic in initial_conditions:
    t, Y = integrate(lv, np.array(ic, dtype=float), dt, T)
    results.append((ic, t, Y))

# --- Plots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

# Time series (first initial condition shown for clarity)
ic0, t0, Y0 = results[0]
ax1.plot(t0, Y0[:, 0], label="Prey N")
ax1.plot(t0, Y0[:, 1], label="Predator P")
ax1.set_xlabel("time")
ax1.set_ylabel("population")
ax1.set_title(f"Time series (IC N,P = {ic0})")
ax1.legend()

# Phase plane: nested closed orbits from each IC
for ic, t, Y in results:
    ax2.plot(Y[:, 0], Y[:, 1], label=f"IC {ic}")
ax2.set_xlabel("Prey N")
ax2.set_ylabel("Predator P")
ax2.set_title("Phase plane: nested closed orbits")
ax2.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.2.1_s5.png")

# --- Checks ---
print(f"Parameters: a={a}, b={b}, c={c}, d={d}")
print(f"dt={dt}, T={T}")

# Coexistence equilibrium (fixed point): N*=d/c, P*=a/b
N_star, P_star = d / c, a / b
print(f"Equilibrium N* = d/c = {N_star}")
print(f"Equilibrium P* = a/b = {P_star}")

# Check 1: out-of-phase oscillation.
# Compare the times at which N and P first peak; a nonzero lag => out of phase.
for ic, t, Y in results:
    iN = int(np.argmax(Y[:, 0]))
    iP = int(np.argmax(Y[:, 1]))
    lag = t[iP] - t[iN]  # predator peak lags prey peak
    # correlation of the two centered signals over the record
    n = Y[:, 0] - Y[:, 0].mean()
    p = Y[:, 1] - Y[:, 1].mean()
    corr = float(np.sum(n * p) / np.sqrt(np.sum(n * n) * np.sum(p * p)))
    print(f"IC {ic}: prey peak t={t[iN]:.2f}, predator peak t={t[iP]:.2f}, "
          f"lag={lag:.2f}, N-P correlation={corr:.3f}")

# Check 2: closed orbits.
# For a closed loop the trajectory returns near its start; report the
# minimum distance back to the start after leaving a neighborhood of it.
for ic, t, Y in results:
    start = Y[0]
    dist = np.linalg.norm(Y - start, axis=1)
    left = np.where(dist > 0.1 * np.linalg.norm(start))[0]
    if left.size:
        after = left[0]
        closure = float(dist[after:].min())
    else:
        closure = float(dist.min())
    print(f"IC {ic}: closest return distance to start = {closure:.4f}")

# One-sentence explanation:
print("Explanation: A near-zero return distance means each trajectory comes back "
      "to its starting point (a closed loop), and the predator peak lagging the "
      "prey peak (negative-to-weak correlation with nonzero lag) means the two "
      "populations oscillate out of phase, together confirming sustained "
      "predator-prey cycles on nested closed orbits.")
