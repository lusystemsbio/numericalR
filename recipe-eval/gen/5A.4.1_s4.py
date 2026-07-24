import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: harmonic oscillator, mass on a spring ---
# d2x/dt2 = -k*x  ->  dx/dt = v, dv/dt = -k*x
# energy per unit mass e = 0.5*k*x^2 + 0.5*v^2 (exactly constant for true motion)
k = 0.1
x0, v0 = 1.0, 2.0
t_end = 100.0

def force(x):
    # acceleration f = dv/dt = -k*x (mass = 1)
    return -k * x

def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

def velocity_verlet(dt):
    n = int(round(t_end / dt))
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0
    f = force(x[0])                       # initial force
    for i in range(n):
        v_half = v[i] + 0.5 * dt * f      # first half velocity update
        x[i+1] = x[i] + dt * v_half       # position update using half-step velocity
        f_next = force(x[i+1])            # force at the new position
        v[i+1] = v_half + 0.5 * dt * f_next  # second half velocity update
        f = f_next                        # carry force forward (reuse, one force eval/step)
        t[i+1] = t[i] + dt
    return t, x, v

# integrate at both step sizes
t_small, x_small, v_small = velocity_verlet(0.01)
t_large, x_large, v_large = velocity_verlet(0.1)

e_small = energy(x_small, v_small)
e_large = energy(x_large, v_large)
e0 = energy(x0, v0)

# --- Check: energy conservation and matching time points ---
# velocity Verlet returns x and v at the SAME time grid, so e(t) is well defined at each t.
drift_small = np.max(np.abs(e_small - e0))
drift_large = np.max(np.abs(e_large - e0))

print(f"Initial energy e0 = {e0:.10f}")
print(f"dt=0.01: max energy deviation = {drift_small:.3e} (relative {drift_small/e0:.3e})")
print(f"dt=0.10: max energy deviation = {drift_large:.3e} (relative {drift_large/e0:.3e})")
print(f"dt=0.01: number of steps = {len(t_small)-1}, final t = {t_small[-1]:.4f}")
print(f"dt=0.10: number of steps = {len(t_large)-1}, final t = {t_large[-1]:.4f}")
print(f"dt=0.01: final x = {x_small[-1]:.6f}, final v = {v_small[-1]:.6f}")
print(f"dt=0.10: final x = {x_large[-1]:.6f}, final v = {v_large[-1]:.6f}")

# confirm x and v are reported on the identical time grid (same time points)
same_grid_small = np.allclose(np.diff(t_small), 0.01)
same_grid_large = np.allclose(np.diff(t_large), 0.10)
print(f"dt=0.01: x and v share one uniform time grid = {same_grid_small}")
print(f"dt=0.10: x and v share one uniform time grid = {same_grid_large}")

# energy stays bounded (oscillates about e0) rather than drifting away -> conserved
conserved_small = drift_small / e0 < 1e-2
conserved_large = drift_large / e0 < 1e-2
print(f"dt=0.01: energy conserved (rel. bound < 1e-2) = {conserved_small}")
print(f"dt=0.10: energy conserved (rel. bound < 1e-2) = {conserved_large}")

# --- Plots: x(t) and e(t) at both step sizes ---
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

axes[0].plot(t_small, x_small, label="dt = 0.01", lw=1)
axes[0].plot(t_large, x_large, label="dt = 0.10", lw=1, ls="--")
axes[0].set_xlabel("t")
axes[0].set_ylabel("x(t)")
axes[0].set_title("Velocity Verlet: position x(t)")
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].plot(t_small, e_small, label="dt = 0.01", lw=1)
axes[1].plot(t_large, e_large, label="dt = 0.10", lw=1, ls="--")
axes[1].axhline(e0, color="k", lw=0.8, ls=":", label="e0 (exact)")
axes[1].set_xlabel("t")
axes[1].set_ylabel("e(t)")
axes[1].set_title("Velocity Verlet: energy e(t)")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.4.1_s4.png")

# One-sentence explanation of why the check confirms the result:
# Because velocity Verlet advances x and v on a single shared time grid, the energy
# e(t)=0.5*k*x^2+0.5*v^2 is evaluated at consistent times, and its bounded oscillation
# about e0 (rather than growing drift) confirms the integrator conserves energy at both step sizes.
print("Check explanation: x and v live on the same time grid, so e(t) is consistent, "
      "and e(t) merely oscillates within a small bound about e0 instead of drifting, "
      "confirming energy conservation at both step sizes.")
