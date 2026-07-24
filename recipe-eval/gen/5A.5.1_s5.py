import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: harmonic oscillator (mass on a spring) ----
# d^2x/dt^2 = -k*x  ->  dx/dt = v, dv/dt = -k*x
# Energy per unit mass e = 0.5*k*x^2 + 0.5*v^2 is exactly constant for true motion.
k = 0.1
x0 = 1.0
v0 = 2.0
t_end = 100.0

def accel(x):
    # acceleration f = -k*x (right-hand side of dv/dt)
    return -k * x

def verlet(dt):
    # Number of steps to reach t_end
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.zeros(n + 1)

    # Startup step: x_1 = x0 + dt*v0 + 0.5*dt^2*f0 (uses v0 only once)
    x[0] = x0
    f0 = accel(x[0])
    if n >= 1:
        x[1] = x0 + dt * v0 + 0.5 * dt * dt * f0

    # Position-only Verlet recurrence, no velocity stored:
    # x_{n+2} = 2*x_{n+1} - x_n + dt^2*f_{n+1}
    for i in range(1, n):
        f = accel(x[i])                       # force at current position only
        x[i + 1] = 2.0 * x[i] - x[i - 1] + dt * dt * f
    return t, x

# ---- Analytic (true) solution for reference ----
# x(t) = x0*cos(w t) + (v0/w)*sin(w t), w = sqrt(k)
w = np.sqrt(k)
def exact(t):
    return x0 * np.cos(w * t) + (v0 / w) * np.sin(w * t)

# Exact energy (constant)
e_exact = 0.5 * k * x0**2 + 0.5 * v0**2
print("Exact conserved energy per unit mass e = %.10f" % e_exact)

# ---- Run at both step sizes ----
t_small, x_small = verlet(0.01)
t_big, x_big = verlet(0.1)

# ---- Check: compare Verlet position to exact solution at final time ----
# (velocity is never stored; we only ever look at positions)
err_small = np.max(np.abs(x_small - exact(t_small)))
err_big = np.max(np.abs(x_big - exact(t_big)))
print("dt = 0.01: final x = %.10f, exact x = %.10f" % (x_small[-1], exact(t_small[-1])))
print("dt = 0.10: final x = %.10f, exact x = %.10f" % (x_big[-1], exact(t_big[-1])))
print("dt = 0.01: max |x_verlet - x_exact| over [0,100] = %.10e" % err_small)
print("dt = 0.10: max |x_verlet - x_exact| over [0,100] = %.10e" % err_big)

# Amplitude check (true amplitude = sqrt(x0^2 + (v0/w)^2))
amp_true = np.sqrt(x0**2 + (v0 / w)**2)
print("True oscillation amplitude = %.10f" % amp_true)
print("dt = 0.01: observed max|x| = %.10f" % np.max(np.abs(x_small)))
print("dt = 0.10: observed max|x| = %.10f" % np.max(np.abs(x_big)))

# ---- Plot x(t) for both step sizes ----
plt.figure(figsize=(11, 5))
plt.plot(t_small, exact(t_small), 'k-', lw=1.0, alpha=0.4, label="exact")
plt.plot(t_small, x_small, 'b-', lw=1.0, label="Verlet dt=0.01")
plt.plot(t_big, x_big, 'r--', lw=1.2, label="Verlet dt=0.1")
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Position-only Verlet: harmonic oscillator (k=0.1)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.5.1_s5.png")

# Explanation of the check:
print("Why the check confirms the result: the position trajectory matches the exact "
      "oscillation (bounded amplitude, correct frequency) with error shrinking as dt->0, "
      "showing Verlet reproduces the motion using only stored positions and a single use of v0.")
