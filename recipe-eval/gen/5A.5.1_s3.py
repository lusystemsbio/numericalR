import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters (harmonic oscillator, mass on a spring) ---
k = 0.1          # spring stiffness
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
T = 100.0        # final time

# acceleration f(x) = -k*x  (d^2x/dt^2 = -k*x)
def f(x):
    return -k * x

# exact energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2 (constant for true motion)
# exact solution of dx/dt=v, dv/dt=-k*x with x(0)=x0, v(0)=v0
w = np.sqrt(k)                       # angular frequency
def x_exact(t):
    return x0 * np.cos(w * t) + (v0 / w) * np.sin(w * t)

# --- Position-only Verlet integrator (stores positions only) ---
def verlet(dt):
    n = int(round(T / dt))           # number of steps
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1)

    # store only positions; velocity is never kept
    x[0] = x0
    # startup step from Taylor expansion: x_1 = x0 + dt*v0 + 0.5*dt^2*f0
    x[1] = x0 + dt * v0 + 0.5 * dt**2 * f(x0)

    # advance position from the two previous positions:
    # x_{n+2} = 2*x_{n+1} - x_n + dt^2 * f_{n+1}
    for i in range(1, n):
        x[i + 1] = 2.0 * x[i] - x[i - 1] + dt**2 * f(x[i])

    return t, x

# run at both step sizes
t_fine, x_fine = verlet(0.01)
t_coarse, x_coarse = verlet(0.1)

# --- Check: compare final position and max abs error against exact solution ---
xe_fine = x_exact(t_fine)
xe_coarse = x_exact(t_coarse)
err_fine = np.max(np.abs(x_fine - xe_fine))
err_coarse = np.max(np.abs(x_coarse - xe_coarse))

# initial energy per unit mass (reference for the true motion)
e0 = 0.5 * k * x0**2 + 0.5 * v0**2

print("Initial energy per unit mass e0 = 0.5*k*x0^2 + 0.5*v0^2 =", e0)
print("Exact final position x_exact(t=100) =", x_exact(T))
print("Verlet dt=0.01: final position x(100) =", x_fine[-1])
print("Verlet dt=0.10: final position x(100) =", x_coarse[-1])
print("Verlet dt=0.01: max |x - x_exact| over [0,100] =", err_fine)
print("Verlet dt=0.10: max |x - x_exact| over [0,100] =", err_coarse)
print("Verlet dt=0.01: amplitude (max|x|) =", np.max(np.abs(x_fine)), "exact amplitude =", np.max(np.abs(xe_fine)))
print("Verlet dt=0.10: amplitude (max|x|) =", np.max(np.abs(x_coarse)), "exact amplitude =", np.max(np.abs(xe_coarse)))

# --- Plot x(t) for Verlet at both step sizes ---
plt.figure(figsize=(11, 5))
plt.plot(t_fine, xe_fine, 'k-', lw=1.0, alpha=0.4, label='exact')
plt.plot(t_fine, x_fine, 'b-', lw=1.0, label='Verlet dt=0.01')
plt.plot(t_coarse, x_coarse, 'r--', lw=1.0, label='Verlet dt=0.1')
plt.xlabel('t')
plt.ylabel('x(t)')
plt.title('Position-only Verlet: harmonic oscillator (k=0.1, x0=1, v0=2)')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.5.1_s3.png")

# Explanation of why the check confirms the result:
# The small max |x - x_exact| at both step sizes shows Verlet tracks the true
# bounded oscillation using only stored positions, confirming it reproduces the
# motion without ever needing the velocity.
print("Why the check confirms it: the small max position error and preserved oscillation amplitude at both step sizes show Verlet reproduces the true bounded motion using stored positions alone (no velocity), which is why it and velocity Verlet underlie most practical MD codes.")
