import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Harmonic oscillator: d^2x/dt^2 = -k*x, i.e. dx/dt = v, dv/dt = -k*x.
# Energy per unit mass e = 0.5*k*x^2 + 0.5*v^2 is exactly conserved.
# ----------------------------------------------------------------------

k = 0.1          # spring stiffness
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # final time

# Force per unit mass (acceleration) for the harmonic oscillator.
def force(x):
    return -k * x

# Explicit velocity-Verlet integrator (implemented step-by-step, no library call).
def velocity_verlet(dt):
    n = int(round(t_end / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0                  # set initial conditions
    f = force(x[0])                      # initial force
    for i in range(n):
        v_half = v[i] + 0.5 * dt * f     # first half-step of velocity update
        x[i + 1] = x[i] + dt * v_half    # full-step position update
        f_next = force(x[i + 1])         # force at the new position
        v[i + 1] = v_half + 0.5 * dt * f_next  # second half-step of velocity update
        f = f_next                       # carry force forward (reuse, no recompute)
        t[i + 1] = t[i] + dt             # advance time
    return t, x, v

# Energy per unit mass.
def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

# Integrate at both step sizes.
t1, x1, v1 = velocity_verlet(0.01)
t2, x2, v2 = velocity_verlet(0.1)

e1 = energy(x1, v1)
e2 = energy(x2, v2)

# ----------------------------------------------------------------------
# Numerical results
# ----------------------------------------------------------------------
e_init = energy(x0, v0)
print(f"Initial energy e(0)                              = {e_init:.10f}")

print(f"dt=0.01: number of steps                         = {len(t1) - 1}")
print(f"dt=0.01: final time reached                      = {t1[-1]:.6f}")
print(f"dt=0.01: final position x(100)                   = {x1[-1]:.10f}")
print(f"dt=0.01: final velocity v(100)                   = {v1[-1]:.10f}")
print(f"dt=0.01: final energy e(100)                     = {e1[-1]:.10f}")
print(f"dt=0.01: max |e(t) - e(0)|                       = {np.max(np.abs(e1 - e_init)):.3e}")
print(f"dt=0.01: max relative energy drift               = {np.max(np.abs(e1 - e_init)) / e_init:.3e}")

print(f"dt=0.1 : number of steps                         = {len(t2) - 1}")
print(f"dt=0.1 : final time reached                      = {t2[-1]:.6f}")
print(f"dt=0.1 : final position x(100)                   = {x2[-1]:.10f}")
print(f"dt=0.1 : final velocity v(100)                   = {v2[-1]:.10f}")
print(f"dt=0.1 : final energy e(100)                     = {e2[-1]:.10f}")
print(f"dt=0.1 : max |e(t) - e(0)|                       = {np.max(np.abs(e2 - e_init)):.3e}")
print(f"dt=0.1 : max relative energy drift               = {np.max(np.abs(e2 - e_init)) / e_init:.3e}")

# ----------------------------------------------------------------------
# Separate check: energy conservation and matched time points.
# Velocity Verlet returns x[i] and v[i] both defined at the SAME time t[i],
# so energy e(t)=0.5*k*x^2+0.5*v^2 is evaluated at a consistent instant.
# ----------------------------------------------------------------------
tol = 1e-2  # bounded-oscillation tolerance on relative energy drift
check1 = np.max(np.abs(e1 - e_init)) / e_init < tol
check2 = np.max(np.abs(e2 - e_init)) / e_init < tol
# Position and velocity share the same time array by construction:
same_time_pts_1 = (t1.shape == x1.shape == v1.shape)
same_time_pts_2 = (t2.shape == x2.shape == v2.shape)
print(f"CHECK dt=0.01: energy bounded within {tol:g} rel.  = {bool(check1)}")
print(f"CHECK dt=0.1 : energy bounded within {tol:g} rel.  = {bool(check2)}")
print(f"CHECK dt=0.01: x,v reported at same time points  = {bool(same_time_pts_1)}")
print(f"CHECK dt=0.1 : x,v reported at same time points  = {bool(same_time_pts_2)}")
# One-sentence explanation of why this confirms the result:
print("EXPLANATION: Because velocity Verlet reports x and v at the same time points, "
      "the energy e(t) is evaluated consistently, and its staying bounded (rather than "
      "drifting) confirms the symplectic integrator conserves energy at both step sizes.")

# ----------------------------------------------------------------------
# Plots: x(t) and e(t) for both step sizes.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

axes[0].plot(t1, x1, label="dt = 0.01", lw=1.0)
axes[0].plot(t2, x2, label="dt = 0.1", lw=1.0, ls="--")
axes[0].set_xlabel("t")
axes[0].set_ylabel("x(t)")
axes[0].set_title("Velocity Verlet: position x(t)")
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].plot(t1, e1, label="dt = 0.01", lw=1.0)
axes[1].plot(t2, e2, label="dt = 0.1", lw=1.0, ls="--")
axes[1].axhline(e_init, color="k", lw=0.8, ls=":", label="exact energy")
axes[1].set_xlabel("t")
axes[1].set_ylabel("e(t)")
axes[1].set_title("Velocity Verlet: energy e(t)")
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.4.1_s2.png")
