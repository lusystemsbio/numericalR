import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Harmonic oscillator: d^2x/dt^2 = -k*x  <=>  dx/dt = v, dv/dt = -k*x
# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2 (constant for the true motion).

k = 0.1
x0, v0 = 1.0, 2.0
t_end = 100.0

def euler(k, x0, v0, dt, t_end):
    n = int(round(t_end / dt))          # number of steps
    t = np.linspace(0.0, n * dt, n + 1) # time grid
    x = np.empty(n + 1)
    v = np.empty(n + 1)
    x[0], v[0] = x0, v0
    for i in range(n):
        f = -k * x[i]                   # force per unit mass, evaluated at current x
        v[i + 1] = v[i] + dt * f        # explicit Euler: update velocity first
        x[i + 1] = x[i] + dt * v[i]     # then update position (uses old velocity)
    e = 0.5 * k * x**2 + 0.5 * v**2     # energy per unit mass along the trajectory
    return t, x, v, e

# Integrate for both step sizes
t_fine, x_fine, v_fine, e_fine = euler(k, x0, v0, 0.01, t_end)
t_coarse, x_coarse, v_coarse, e_coarse = euler(k, x0, v0, 0.1, t_end)

# Initial (true, constant) energy
e_true = 0.5 * k * x0**2 + 0.5 * v0**2
print(f"Initial/true energy e0                : {e_true:.6f}")

# Report energies and drift at the two step sizes
print(f"dt=0.01: final energy e(100)          : {e_fine[-1]:.6f}")
print(f"dt=0.01: energy drift e(100)-e0       : {e_fine[-1] - e_true:.6f}")
print(f"dt=0.01: relative drift               : {(e_fine[-1] - e_true) / e_true:.6f}")
print(f"dt=0.01: max |x| (amplitude)          : {np.max(np.abs(x_fine)):.6f}")

print(f"dt=0.1 : final energy e(100)          : {e_coarse[-1]:.6f}")
print(f"dt=0.1 : energy drift e(100)-e0       : {e_coarse[-1] - e_true:.6f}")
print(f"dt=0.1 : relative drift               : {(e_coarse[-1] - e_true) / e_true:.6f}")
print(f"dt=0.1 : max |x| (amplitude)          : {np.max(np.abs(x_coarse)):.6f}")

# Confirm monotonic upward energy drift (steady growth)
diffs_coarse = np.diff(e_coarse)
diffs_fine = np.diff(e_fine)
print(f"dt=0.1 : energy strictly increasing?  : {bool(np.all(diffs_coarse > 0))}")
print(f"dt=0.01: energy strictly increasing?  : {bool(np.all(diffs_fine > 0))}")

# Plots: x(t) and e(t) for both step sizes
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0, 0].plot(t_fine, x_fine, lw=0.8)
axes[0, 0].set_title("x(t), dt = 0.01")
axes[0, 0].set_xlabel("t"); axes[0, 0].set_ylabel("x")

axes[0, 1].plot(t_coarse, x_coarse, lw=0.8, color="C1")
axes[0, 1].set_title("x(t), dt = 0.1 (amplitude grows)")
axes[0, 1].set_xlabel("t"); axes[0, 1].set_ylabel("x")

axes[1, 0].plot(t_fine, e_fine, lw=0.8)
axes[1, 0].axhline(e_true, color="k", ls="--", lw=0.8, label="true e")
axes[1, 0].set_title("e(t), dt = 0.01 (slow upward drift)")
axes[1, 0].set_xlabel("t"); axes[1, 0].set_ylabel("e"); axes[1, 0].legend()

axes[1, 1].plot(t_coarse, e_coarse, lw=0.8, color="C1")
axes[1, 1].axhline(e_true, color="k", ls="--", lw=0.8, label="true e")
axes[1, 1].set_title("e(t), dt = 0.1 (energy grows)")
axes[1, 1].set_xlabel("t"); axes[1, 1].set_ylabel("e"); axes[1, 1].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.2.1_s1.png")

# Explanation of why the check confirms non-conservation:
print("Explanation: Because the true energy is exactly constant, any steady rise in "
      "the computed e(t)-larger at dt=0.1, smaller but still positive at dt=0.01-"
      "shows the explicit Euler scheme injects spurious energy and does not conserve it.")
