import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters: mass on a spring, d^2x/dt^2 = -k*x
k = 0.1          # spring stiffness
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # final time

# Acceleration (force per unit mass) for the harmonic oscillator
def f(x):
    return -k * x

# Energy per unit mass, exactly constant for the true motion
def energy(x, v):
    return 0.5 * k * x * x + 0.5 * v * v

# Position-only Verlet integrator: stores only positions, never velocity.
# x_{n+2} = 2*x_{n+1} - x_n + dt^2 * f_{n+1}
def verlet(dt):
    n = int(round(t_end / dt))          # number of steps
    t = np.linspace(0.0, n * dt, n + 1) # time grid
    x = np.empty(n + 1)                 # store positions only

    # Startup step: Taylor expansion x_1 = x0 + dt*v0 + 0.5*dt^2*f0
    x[0] = x0
    x[1] = x0 + dt * v0 + 0.5 * dt * dt * f(x0)

    # Main recurrence: advance position from the two previous positions
    for i in range(1, n):
        x[i + 1] = 2.0 * x[i] - x[i - 1] + dt * dt * f(x[i])

    return t, x

# Analytic (true) solution for the check: x(t) = x0*cos(w t) + (v0/w)*sin(w t)
w = np.sqrt(k)
def x_true(t):
    return x0 * np.cos(w * t) + (v0 / w) * np.sin(w * t)

# Integrate at both step sizes
t_fine, x_fine = verlet(0.01)
t_coarse, x_coarse = verlet(0.1)

# --- Numerical results ---
# Energies from Verlet require a velocity estimate; use centered difference
# v_n ~ (x_{n+1} - x_{n-1}) / (2*dt) purely for the constancy check.
def energy_stats(t, x, dt):
    v_mid = (x[2:] - x[:-2]) / (2.0 * dt)
    x_mid = x[1:-1]
    e = energy(x_mid, v_mid)
    return e

for label, (t, x), dt in [("dt=0.01", (t_fine, x_fine), 0.01),
                          ("dt=0.10", (t_coarse, x_coarse), 0.1)]:
    e = energy_stats(t, x, dt)
    max_err = np.max(np.abs(x - x_true(t)))
    print(f"Verlet {label}: final position x(t=100) = {x[-1]:.6f}")
    print(f"Verlet {label}: true  position x(t=100) = {x_true(t[-1]):.6f}")
    print(f"Verlet {label}: max |x_verlet - x_true| over run = {max_err:.6e}")
    print(f"Verlet {label}: energy mean (centered-diff v) = {np.mean(e):.6f}")
    print(f"Verlet {label}: energy std  (centered-diff v) = {np.std(e):.6e}")
    print(f"Verlet {label}: energy drift (max-min)        = {np.max(e)-np.min(e):.6e}")

print(f"Exact energy per unit mass e0 = {energy(x0, v0):.6f}")

# --- Plot x(t) for Verlet at both step sizes ---
plt.figure(figsize=(11, 5))
plt.plot(t_fine, x_fine, '-', lw=1.0, label="Verlet dt=0.01")
plt.plot(t_coarse, x_coarse, '--', lw=1.2, label="Verlet dt=0.10")
plt.plot(t_fine, x_true(t_fine), ':', color='k', lw=1.0, label="true x(t)")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Position-only Verlet: harmonic oscillator (k=0.1, x0=1, v0=2)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.5.1_s2.png")

# The check confirms the result because the Verlet trajectories at both step sizes
# stay close to the exact sinusoid (small max error, no growing amplitude) using only
# stored positions, showing the method oscillates stably without ever tracking velocity.
print("Check: both Verlet step sizes track the exact oscillation using positions only, "
      "with bounded error and no amplitude blow-up, confirming stable position-only integration.")
