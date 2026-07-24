import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
k = 0.1           # spring stiffness
x0 = 1.0          # initial position
v0 = 2.0          # initial velocity
t_end = 100.0     # final time

# Acceleration (force per unit mass) for the harmonic oscillator
def f(x):
    return -k * x

# Energy per unit mass: exactly constant for the true motion
def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

# --- Explicit (forward) Euler integrator, for comparison ---
def euler(dt):
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1)
    v = np.empty(n + 1)
    x[0], v[0] = x0, v0
    for i in range(n):
        # update both variables using values at the current (same) time
        x[i + 1] = x[i] + dt * v[i]
        v[i + 1] = v[i] + dt * f(x[i])
    return t, x, v

# --- Leapfrog integrator with velocity carried at half-integer steps ---
def leapfrog(dt):
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1)
    v = np.empty(n + 1)   # velocity synced back to integer steps (for reporting/energy)
    x[0], v[0] = x0, v0

    # startup half-step: kick the velocity forward by half a step
    v_half = v0 + 0.5 * dt * f(x0)

    for i in range(n):
        # drift: advance position a full step using the half-step velocity
        x[i + 1] = x[i] + dt * v_half
        # kick: advance the half-step velocity a full step using the new force
        v_next_half = v_half + dt * f(x[i + 1])
        # synchronized velocity at the integer step is the average of the two half-steps
        v[i + 1] = 0.5 * (v_half + v_next_half)
        # advance the half-step velocity for the next iteration
        v_half = v_next_half
    return t, x, v

# --- Run all cases ---
dts = [0.01, 0.1]
results = {}
for dt in dts:
    t_lf, x_lf, v_lf = leapfrog(dt)
    t_eu, x_eu, v_eu = euler(dt)
    results[dt] = dict(t_lf=t_lf, x_lf=x_lf, v_lf=v_lf,
                       t_eu=t_eu, x_eu=x_eu, v_eu=v_eu)

# --- Numerical reporting ---
e_true = energy(x0, v0)
print(f"Exact conserved energy e0 = {e_true:.8f}")

for dt in dts:
    r = results[dt]
    e_lf = energy(r["x_lf"], r["v_lf"])
    e_eu = energy(r["x_eu"], r["v_eu"])
    print(f"--- dt = {dt} ---")
    print(f"Leapfrog: final x = {r['x_lf'][-1]:.6f}, final v = {r['v_lf'][-1]:.6f}")
    print(f"Leapfrog: energy min = {e_lf.min():.8f}, max = {e_lf.max():.8f}")
    print(f"Leapfrog: energy drift (final - initial) = {e_lf[-1] - e_lf[0]:.8e}")
    print(f"Leapfrog: max |e - e0| over run          = {np.max(np.abs(e_lf - e_true)):.8e}")
    print(f"Euler:    final x = {r['x_eu'][-1]:.6f}, final v = {r['v_eu'][-1]:.6f}")
    print(f"Euler:    energy min = {e_eu.min():.8f}, max = {e_eu.max():.8f}")
    print(f"Euler:    energy drift (final - initial) = {e_eu[-1] - e_eu[0]:.8e}")

# The check confirms the result because leapfrog's energy error stays bounded (a
# small bounded ripple) instead of growing without limit like Euler's: a bounded
# error with no secular trend is exactly the signature of a symplectic, time-
# reversible integrator conserving energy over long times.

# --- Plots ---
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

for j, dt in enumerate(dts):
    r = results[dt]
    e_lf = energy(r["x_lf"], r["v_lf"])
    e_eu = energy(r["x_eu"], r["v_eu"])

    ax_x = axes[0, j]
    ax_x.plot(r["t_lf"], r["x_lf"], label="leapfrog", lw=1.0)
    ax_x.plot(r["t_eu"], r["x_eu"], label="Euler", lw=0.8, alpha=0.7)
    ax_x.set_title(f"x(t), dt = {dt}")
    ax_x.set_xlabel("t"); ax_x.set_ylabel("x")
    ax_x.legend()

    ax_e = axes[1, j]
    ax_e.plot(r["t_lf"], e_lf, label="leapfrog", lw=1.0)
    ax_e.plot(r["t_eu"], e_eu, label="Euler", lw=0.8, alpha=0.7)
    ax_e.axhline(e_true, color="k", ls="--", lw=0.8, label="exact e0")
    ax_e.set_title(f"energy e(t), dt = {dt}")
    ax_e.set_xlabel("t"); ax_e.set_ylabel("e")
    ax_e.legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.3.1_s5.png", dpi=120)
