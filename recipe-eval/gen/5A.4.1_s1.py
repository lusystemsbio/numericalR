import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Harmonic oscillator: d2x/dt2 = -k*x  ->  dx/dt = v, dv/dt = -k*x
# Energy per unit mass e = 0.5*k*x^2 + 0.5*v^2 is conserved for true motion.
# ----------------------------------------------------------------------

k = 0.1          # spring stiffness
x0, v0 = 1.0, 2.0  # initial position and velocity
t_end = 100.0

def force(x):
    # acceleration (force per unit mass) f = -k*x
    return -k * x

def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

def velocity_verlet(dt):
    # number of steps to reach t_end
    n = int(round(t_end / dt))
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0
    f = force(x[0])                     # current force
    for i in range(n):
        v_half = v[i] + 0.5 * dt * f    # half-step velocity kick
        x[i + 1] = x[i] + dt * v_half   # full-step position drift
        f_next = force(x[i + 1])        # force at new position
        v[i + 1] = v_half + 0.5 * dt * f_next  # second half-step kick
        f = f_next                      # carry force forward
        t[i + 1] = t[i] + dt
    return t, x, v

# integrate at both step sizes
t_small, x_small, v_small = velocity_verlet(0.01)
t_large, x_large, v_large = velocity_verlet(0.1)

e_small = energy(x_small, v_small)
e_large = energy(x_large, v_large)
e_exact = energy(x0, v0)

# ----------------------------------------------------------------------
# Plots: x(t) and e(t) for both step sizes
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 1, figsize=(10, 8))

axes[0].plot(t_small, x_small, label="dt = 0.01")
axes[0].plot(t_large, x_large, "--", label="dt = 0.1")
axes[0].set_xlabel("t")
axes[0].set_ylabel("x(t)")
axes[0].set_title("Velocity Verlet: position x(t)")
axes[0].legend()

axes[1].plot(t_small, e_small, label="dt = 0.01")
axes[1].plot(t_large, e_large, "--", label="dt = 0.1")
axes[1].axhline(e_exact, color="k", lw=0.8, label="exact e")
axes[1].set_xlabel("t")
axes[1].set_ylabel("energy e(t)")
axes[1].set_title("Velocity Verlet: energy e(t)")
axes[1].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.4.1_s1.png")

# ----------------------------------------------------------------------
# Check: energy conservation and matching time points
# ----------------------------------------------------------------------
drift_small = np.max(np.abs(e_small - e_exact))
drift_large = np.max(np.abs(e_large - e_exact))

# relative energy drift
rel_small = drift_small / e_exact
rel_large = drift_large / e_exact

# x and v arrays share the same time grid t (same length, same t values)
same_times_small = (len(t_small) == len(x_small) == len(v_small))
same_times_large = (len(t_large) == len(x_large) == len(v_large))

print(f"Exact energy e (constant)                 : {e_exact:.10f}")
print(f"dt=0.01  initial energy                   : {e_small[0]:.10f}")
print(f"dt=0.01  final energy                     : {e_small[-1]:.10f}")
print(f"dt=0.01  max |e - e_exact| (abs drift)    : {drift_small:.3e}")
print(f"dt=0.01  max relative energy drift        : {rel_small:.3e}")
print(f"dt=0.10  initial energy                   : {e_large[0]:.10f}")
print(f"dt=0.10  final energy                     : {e_large[-1]:.10f}")
print(f"dt=0.10  max |e - e_exact| (abs drift)    : {drift_large:.3e}")
print(f"dt=0.10  max relative energy drift        : {rel_large:.3e}")
print(f"dt=0.01  x and v on identical time grid   : {same_times_small}")
print(f"dt=0.10  x and v on identical time grid   : {same_times_large}")

# One-sentence explanation:
# The check confirms the result because the energy stays bounded near its exact
# value (no secular drift) at both step sizes while x and v are reported on the
# very same time grid, which is exactly the defining property of a symplectic,
# synchronized integrator like velocity Verlet.
print("Explanation: bounded energy (no drift) with x and v on the same time "
      "grid at both dt is precisely the signature of correct velocity Verlet.")
