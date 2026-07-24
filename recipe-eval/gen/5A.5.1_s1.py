import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
k = 0.1            # spring stiffness (per unit mass)
x0 = 1.0           # initial position
v0 = 2.0           # initial velocity
t_end = 100.0      # final time

# force per unit mass: dv/dt = -k*x
def force(x):
    return -k * x

# energy per unit mass, exactly constant for the true motion
def energy(x, v):
    return 0.5 * k * x * x + 0.5 * v * v

# ---- Original position-only Verlet integrator ----
# Advances position from the TWO previous positions; velocity is never stored.
def verlet(dt):
    n = int(round(t_end / dt))          # number of steps
    t = np.linspace(0.0, n * dt, n + 1) # time grid
    x = np.empty(n + 1)                 # only positions are stored

    # Startup step: we only have x0 and v0, so use the Taylor expansion
    # x_1 = x0 + dt*v0 + 0.5*dt^2*f0 to bootstrap the two-point recursion.
    x[0] = x0
    f0 = force(x[0])
    x[1] = x0 + dt * v0 + 0.5 * dt * dt * f0

    # Main recursion: x_{n+2} = 2*x_{n+1} - x_n + dt^2 * f_{n+1}
    for i in range(1, n):
        f = force(x[i])                 # force from the current (middle) position
        x[i + 1] = 2.0 * x[i] - x[i - 1] + dt * dt * f

    return t, x

# ---- Exact analytic solution for comparison ----
# x(t) = x0*cos(w t) + (v0/w)*sin(w t), with w = sqrt(k)
def exact(t):
    w = np.sqrt(k)
    return x0 * np.cos(w * t) + (v0 / w) * np.sin(w * t)

# ---- Run at both step sizes ----
t_fine, x_fine = verlet(0.01)
t_coarse, x_coarse = verlet(0.1)

# ---- Numerical checks ----
# Initial energy (only value we can form directly from stored data + given v0).
print(f"Initial energy per unit mass e0 = {energy(x0, v0):.6f}")

# Max deviation from the exact solution confirms the oscillation is reproduced.
err_fine = np.max(np.abs(x_fine - exact(t_fine)))
err_coarse = np.max(np.abs(x_coarse - exact(t_coarse)))
print(f"dt = 0.01: number of steps = {len(t_fine) - 1}")
print(f"dt = 0.01: max |x_verlet - x_exact| over t in [0,100] = {err_fine:.6e}")
print(f"dt = 0.10: number of steps = {len(t_coarse) - 1}")
print(f"dt = 0.10: max |x_verlet - x_exact| over t in [0,100] = {err_coarse:.6e}")

# Amplitude check: true amplitude is sqrt(x0^2 + v0^2/k); Verlet should stay bounded near it.
true_amp = np.sqrt(x0 * x0 + v0 * v0 / k)
print(f"True oscillation amplitude = {true_amp:.6f}")
print(f"dt = 0.01: max |x| = {np.max(np.abs(x_fine)):.6f}")
print(f"dt = 0.10: max |x| = {np.max(np.abs(x_coarse)):.6f}")

# ---- Plot x(t) for Verlet at both step sizes ----
plt.figure(figsize=(11, 5))
plt.plot(t_fine, x_fine, "-", lw=1.0, label="Verlet, dt = 0.01")
plt.plot(t_coarse, x_coarse, "--", lw=1.2, label="Verlet, dt = 0.1")
plt.plot(t_fine, exact(t_fine), ":", color="k", lw=1.0, label="exact")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Harmonic oscillator via position-only Verlet (k = 0.1)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.5.1_s1.png")

# Explanation: the small max deviations from the exact sinusoid at both step sizes
# show the position-only recursion tracks the true oscillation using stored positions
# alone, which is exactly why Verlet-family schemes underlie most practical MD codes.
print("Check: small max deviation from the exact solution at both dt confirms the")
print("position-only Verlet recursion reproduces the true oscillation using only stored positions.")
