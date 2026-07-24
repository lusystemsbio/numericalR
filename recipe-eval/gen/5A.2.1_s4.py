import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
k = 0.1          # spring stiffness (per unit mass)
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # final time

def energy(x, v):
    # energy per unit mass: 0.5*k*x^2 + 0.5*v^2, exactly constant for true motion
    return 0.5 * k * x**2 + 0.5 * v**2

def euler_oscillator(dt):
    """Explicit (forward) Euler for dx/dt = v, dv/dt = -k*x.
    Update velocity first, then position, both using OLD values."""
    n = int(round(t_end / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0                  # set initial condition
    for i in range(n):
        f = -k * x[i]                    # force per unit mass at current position
        v[i+1] = v[i] + dt * f           # v_next = v + dt*f  (velocity update)
        x[i+1] = x[i] + dt * v[i]        # x_next = x + dt*v  (position update, old v)
        t[i+1] = t[i] + dt               # advance time
    return t, x, v

# ---- Run both step sizes ----
results = {}
for dt in (0.01, 0.1):
    t, x, v = euler_oscillator(dt)
    e = energy(x, v)
    results[dt] = (t, x, v, e)

# ---- Print numerical diagnostics ----
e_true = energy(x0, v0)
print(f"Exact conserved energy (per unit mass) e0 = {e_true:.6f}")
for dt in (0.01, 0.1):
    t, x, v, e = results[dt]
    print(f"--- dt = {dt} ---")
    print(f"dt = {dt}: energy at t=0        = {e[0]:.6f}")
    print(f"dt = {dt}: energy at t={t_end:g}       = {e[-1]:.6f}")
    print(f"dt = {dt}: energy drift (end-start) = {e[-1]-e[0]:.6f}")
    print(f"dt = {dt}: relative energy growth   = {(e[-1]-e[0])/e[0]*100:.4f} %")
    print(f"dt = {dt}: max |x| (amplitude)      = {np.max(np.abs(x)):.6f}")

# ---- Plots ----
fig, axes = plt.subplots(2, 2, figsize=(12, 8))
for col, dt in enumerate((0.01, 0.1)):
    t, x, v, e = results[dt]
    axes[0, col].plot(t, x, lw=0.8)
    axes[0, col].set_title(f"Position x(t), dt = {dt}")
    axes[0, col].set_xlabel("t"); axes[0, col].set_ylabel("x")
    axes[1, col].plot(t, e, color="C3", lw=0.8)
    axes[1, col].axhline(e_true, color="k", ls="--", lw=0.8, label="true energy")
    axes[1, col].set_title(f"Energy e(t), dt = {dt}")
    axes[1, col].set_xlabel("t"); axes[1, col].set_ylabel("e")
    axes[1, col].legend()
fig.suptitle("Forward Euler on the harmonic oscillator (k=0.1)")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.2.1_s4.png")

# ---- Explanation ----
# The true motion keeps e exactly constant, so any monotonic upward drift in the
# computed e(t) (large and obvious at dt=0.1, small but nonzero at dt=0.01) is
# purely a numerical artifact, confirming that forward Euler does not conserve energy.
print("Explanation: because the exact energy is constant, the steady upward drift in the")
print("computed e(t) at dt=0.1 (and the slow drift at dt=0.01) can only be numerical error,")
print("confirming that forward Euler fails to conserve energy.")
