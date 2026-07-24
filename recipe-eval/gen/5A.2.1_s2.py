import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters -------------------------------------------------------
k = 0.1          # spring stiffness (per unit mass)
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # integrate to this time

def energy(x, v):
    # energy per unit mass: kinetic + potential
    return 0.5 * k * x**2 + 0.5 * v**2

def euler_oscillator(dt):
    # explicit forward Euler for dx/dt = v, dv/dt = -k*x
    n = int(round(t_end / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0                  # initial conditions
    for i in range(n):
        f = -k * x[i]                    # acceleration at current state
        v[i+1] = v[i] + dt * f           # update velocity first
        x[i+1] = x[i] + dt * v[i]        # then update position (uses old v)
        t[i+1] = t[i] + dt               # advance time
    e = energy(x, v)                     # energy at every step
    return t, x, v, e

# --- Run both step sizes ----------------------------------------------------
results = {}
for dt in (0.01, 0.1):
    t, x, v, e = euler_oscillator(dt)
    results[dt] = (t, x, v, e)
    e0, ef = e[0], e[-1]
    print(f"dt = {dt}: initial energy e(0) = {e0:.6f}")
    print(f"dt = {dt}: final energy   e(t_end) = {ef:.6f}")
    print(f"dt = {dt}: energy change  e(t_end) - e(0) = {ef - e0:.6f}")
    print(f"dt = {dt}: relative energy growth = {(ef - e0)/e0*100:.4f} %")
    print(f"dt = {dt}: initial amplitude |x(0)| = {abs(x[0]):.6f}")
    print(f"dt = {dt}: max |x| over run = {np.max(np.abs(x)):.6f}")

# --- Plots ------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

for col, dt in enumerate((0.01, 0.1)):
    t, x, v, e = results[dt]
    axes[0, col].plot(t, x, lw=0.8)
    axes[0, col].set_title(f"x(t), Euler, dt = {dt}")
    axes[0, col].set_xlabel("t")
    axes[0, col].set_ylabel("x")

    axes[1, col].plot(t, e, color="tab:red", lw=0.8)
    axes[1, col].set_title(f"energy e(t), Euler, dt = {dt}")
    axes[1, col].set_xlabel("t")
    axes[1, col].set_ylabel("e")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.2.1_s2.png")

# --- Conservation check -----------------------------------------------------
# The check compares e(t_end) with e(0): the true motion has e exactly
# constant, so any steady positive drift in e (large at dt=0.1, small but
# nonzero at dt=0.01) shows explicit Euler injects energy and does NOT conserve it.
print("Conclusion: e(t_end) > e(0) for both step sizes, growing with dt,")
print("confirming that explicit Euler does not conserve energy.")
