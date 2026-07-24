import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
k = 0.1          # spring stiffness
x0, v0 = 1.0, 2.0  # initial position and velocity
t_end = 100.0    # final time

# Acceleration (force per unit mass) for the harmonic oscillator: dv/dt = -k*x
def accel(x):
    return -k * x

# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2 (exactly constant for true motion)
def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

# Explicit (forward) Euler integrator, written out step by step.
def euler_integrate(dt):
    n = int(round(t_end / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0                  # set initial conditions
    for i in range(n):
        f = accel(x[i])                 # evaluate force at current position
        v[i + 1] = v[i] + dt * f        # update velocity first: v_next = v + dt*f
        x[i + 1] = x[i] + dt * v[i]     # then position: x_next = x + dt*v (old v)
        t[i + 1] = t[i] + dt            # advance time
    e = energy(x, v)                    # energy along the trajectory
    return t, x, v, e

# Run for both step sizes
t_small, x_small, v_small, e_small = euler_integrate(0.01)
t_large, x_large, v_large, e_large = euler_integrate(0.1)

e_init = energy(x0, v0)

# --- Report numerical results ---
print(f"Initial energy e(0)                      = {e_init:.6f}")
print(f"dt=0.01: final energy e(100)             = {e_small[-1]:.6f}")
print(f"dt=0.01: energy drift e(100)-e(0)        = {e_small[-1] - e_init:.6f}")
print(f"dt=0.01: relative energy drift           = {(e_small[-1] - e_init) / e_init:.6f}")
print(f"dt=0.01: max |x|                          = {np.max(np.abs(x_small)):.6f}")
print(f"dt=0.1:  final energy e(100)             = {e_large[-1]:.6f}")
print(f"dt=0.1:  energy drift e(100)-e(0)        = {e_large[-1] - e_init:.6f}")
print(f"dt=0.1:  relative energy drift           = {(e_large[-1] - e_init) / e_init:.6f}")
print(f"dt=0.1:  max |x|                          = {np.max(np.abs(x_large)):.6f}")

# --- Plots ---
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0, 0].plot(t_small, x_small, lw=0.8)
axes[0, 0].set_title("x(t), dt = 0.01")
axes[0, 0].set_xlabel("t"); axes[0, 0].set_ylabel("x")

axes[0, 1].plot(t_large, x_large, lw=0.8, color="C1")
axes[0, 1].set_title("x(t), dt = 0.1")
axes[0, 1].set_xlabel("t"); axes[0, 1].set_ylabel("x")

axes[1, 0].plot(t_small, e_small, lw=0.8)
axes[1, 0].axhline(e_init, color="k", ls="--", lw=0.8, label="e(0)")
axes[1, 0].set_title("energy e(t), dt = 0.01")
axes[1, 0].set_xlabel("t"); axes[1, 0].set_ylabel("e"); axes[1, 0].legend()

axes[1, 1].plot(t_large, e_large, lw=0.8, color="C1")
axes[1, 1].axhline(e_init, color="k", ls="--", lw=0.8, label="e(0)")
axes[1, 1].set_title("energy e(t), dt = 0.1")
axes[1, 1].set_xlabel("t"); axes[1, 1].set_ylabel("e"); axes[1, 1].legend()

fig.suptitle("Forward Euler on the harmonic oscillator: energy is not conserved")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.2.1_s3.png")

# One-sentence explanation:
# Because the true energy is exactly constant, any steady upward growth in the
# computed e(t) can only come from the integrator, so the rising energy at dt=0.1
# (and slow drift at dt=0.01) confirms that forward Euler fails to conserve energy.
print("Explanation: since the exact energy is constant, the observed steady upward "
      "growth/drift in e(t) can only originate from the numerical scheme, confirming "
      "that forward Euler does not conserve energy.")
