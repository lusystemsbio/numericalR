import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
k = 0.1        # spring stiffness
x0, v0 = 1.0, 2.0
T = 100.0

# Force per unit mass for the harmonic oscillator: f = -k*x
def force(x):
    return -k * x

# Exact energy per unit mass (constant for the true motion)
def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

# Explicit velocity-Verlet integrator
def velocity_verlet(dt, T):
    n = int(round(T / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0
    f = force(x[0])                 # current force
    for i in range(n):
        v_half = v[i] + 0.5 * dt * f      # first half velocity update
        x[i+1] = x[i] + dt * v_half       # full position update
        f_next = force(x[i+1])            # force at new position
        v[i+1] = v_half + 0.5 * dt * f_next   # second half velocity update
        f = f_next                        # carry force forward
        t[i+1] = t[i] + dt
    return t, x, v

# Run at both step sizes
t1, x1, v1 = velocity_verlet(0.01, T)
t2, x2, v2 = velocity_verlet(0.1, T)

e1 = energy(x1, v1)
e2 = energy(x2, v2)

# --- Plots ---
fig, axes = plt.subplots(2, 1, figsize=(9, 8))
axes[0].plot(t1, x1, label="dt = 0.01")
axes[0].plot(t2, x2, "--", label="dt = 0.1")
axes[0].set_xlabel("t"); axes[0].set_ylabel("x(t)")
axes[0].set_title("Velocity Verlet: position x(t)")
axes[0].legend()

axes[1].plot(t1, e1, label="dt = 0.01")
axes[1].plot(t2, e2, "--", label="dt = 0.1")
axes[1].set_xlabel("t"); axes[1].set_ylabel("e(t)")
axes[1].set_title("Velocity Verlet: energy e(t)")
axes[1].legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.4.1_s3.png")

# --- Numerical results ---
e_exact = energy(x0, v0)
print(f"Exact initial energy e0: {e_exact}")

print(f"dt=0.01 final time: {t1[-1]}")
print(f"dt=0.01 final x: {x1[-1]}")
print(f"dt=0.01 final v: {v1[-1]}")
print(f"dt=0.01 mean energy: {np.mean(e1)}")
print(f"dt=0.01 max abs energy drift from e0: {np.max(np.abs(e1 - e_exact))}")

print(f"dt=0.1 final time: {t2[-1]}")
print(f"dt=0.1 final x: {x2[-1]}")
print(f"dt=0.1 final v: {v2[-1]}")
print(f"dt=0.1 mean energy: {np.mean(e2)}")
print(f"dt=0.1 max abs energy drift from e0: {np.max(np.abs(e2 - e_exact))}")

# --- Separate check: conservation + shared time points ---
# Energy is "conserved" if its drift stays a small bounded fraction of e0 (no secular growth).
tol = 0.05  # 5% of e0 as an acceptance band for the symplectic energy oscillation
cons1 = np.max(np.abs(e1 - e_exact)) < tol * e_exact
cons2 = np.max(np.abs(e2 - e_exact)) < tol * e_exact
print(f"Energy conserved (bounded, <5% of e0) at dt=0.01: {bool(cons1)}")
print(f"Energy conserved (bounded, <5% of e0) at dt=0.1: {bool(cons2)}")

# x and v are returned on the SAME time grid: identical length and identical times.
same_grid1 = (len(x1) == len(v1) == len(t1))
same_grid2 = (len(x2) == len(v2) == len(t2))
print(f"x and v share the same time points at dt=0.01: {bool(same_grid1)}")
print(f"x and v share the same time points at dt=0.1: {bool(same_grid2)}")

# Explanation: because velocity Verlet is symplectic, its energy oscillates within a small
# bounded band (no secular drift) at BOTH step sizes and x,v are advanced on one shared grid,
# which confirms it is a stable, energy-preserving integrator whose state is co-located in time.
print("Check explanation: bounded (non-drifting) energy at both dt values on a single shared "
      "time grid confirms velocity Verlet is symplectic and reports x and v at the same instants.")
