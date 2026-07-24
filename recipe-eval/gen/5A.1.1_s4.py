import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Test system: harmonic oscillator (mass on a spring)
#   d^2x/dt^2 = -k*x   ->   dx/dt = v ,  dv/dt = -k*x
# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2 (constant for true motion)
# We use direct force evaluation f(x) = -k*x; integrators come later.
# A correct simulation MUST preserve: (1) a pure sinusoidal, bounded
# amplitude, and (2) constant energy e.  Amplitude growth or energy
# drift is therefore purely a numerical artifact.
# ---------------------------------------------------------------

# Force law by direct evaluation (no library ODE routine)
def f(x, k):
    return -k * x  # spring restoring force per unit mass

# Parameters and initial conditions
k = 0.1
x0 = 1.0
v0 = 2.0

print("k (stiffness) =", k)
print("x0 (initial displacement) =", x0)
print("v0 (initial velocity) =", v0)

# ---------------------------------------------------------------
# Exact analytic solution, derived explicitly (not via a solver call).
# General solution of x'' = -k x is x(t) = A*cos(w t) + B*sin(w t),
# with angular frequency w = sqrt(k).
#   x(0) = A         -> A = x0
#   x'(0) = w*B      -> B = v0 / w
# ---------------------------------------------------------------
w = np.sqrt(k)                 # angular frequency
period = 2.0 * np.pi / w       # oscillation period
A = x0                         # cosine coefficient from x(0)
B = v0 / w                     # sine coefficient from v(0)
amplitude = np.sqrt(A**2 + B**2)  # sinusoid amplitude

print("angular frequency w = sqrt(k) =", w)
print("period T = 2*pi/w =", period)
print("cosine coefficient A =", A)
print("sine coefficient B =", B)
print("amplitude sqrt(A^2+B^2) =", amplitude)

# Time grid spanning several periods
t = np.linspace(0.0, 4.0 * period, 2000)

# Build exact motion explicitly, term by term
x_exact = A * np.cos(w * t) + B * np.sin(w * t)          # displacement
v_exact = -A * w * np.sin(w * t) + B * w * np.cos(w * t)  # velocity = dx/dt

# Exact energy per unit mass at every time
e_exact = 0.5 * k * x_exact**2 + 0.5 * v_exact**2

# Initial (reference) energy
e0 = 0.5 * k * x0**2 + 0.5 * v0**2
print("initial energy e0 = 0.5*k*x0^2 + 0.5*v0^2 =", e0)

# ---------------------------------------------------------------
# Check: the exact motion is a pure sinusoid at CONSTANT energy.
# We verify amplitude does not grow (max |x| equals the analytic amplitude)
# and energy does not drift (max variation is at the round-off level).
# ---------------------------------------------------------------
amp_measured = np.max(np.abs(x_exact))          # observed peak displacement
energy_drift = np.max(np.abs(e_exact - e0))     # max deviation from e0
energy_rel_drift = energy_drift / e0            # relative drift

print("measured peak |x| over the run =", amp_measured)
print("amplitude match (analytic - measured) =", amplitude - amp_measured)
print("max absolute energy drift =", energy_drift)
print("max relative energy drift =", energy_rel_drift)

# Explanation (one sentence): because the measured peak displacement equals
# the analytic amplitude and the energy stays flat to round-off, the exact
# motion is confirmed as a pure constant-energy sinusoid, so any amplitude
# growth or energy drift seen later must come from the integrator, not the physics.
print("Check: measured amplitude matches analytic amplitude and energy is "
      "flat to round-off, confirming a pure constant-energy sinusoid; hence "
      "any later amplitude growth or energy drift is a numerical artifact.")

# ---------------------------------------------------------------
# Reference plot: exact motion x(t) and its constant energy e(t)
# ---------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

ax1.plot(t, x_exact, color="tab:blue", label="exact x(t)")
ax1.axhline(amplitude, color="gray", ls="--", lw=0.8, label="+/- amplitude")
ax1.axhline(-amplitude, color="gray", ls="--", lw=0.8)
ax1.set_ylabel("displacement x")
ax1.set_title("Harmonic oscillator (k=%.3g): exact motion and energy reference" % k)
ax1.legend(loc="upper right")
ax1.grid(True, alpha=0.3)

ax2.plot(t, e_exact, color="tab:red", label="exact energy e(t)")
ax2.axhline(e0, color="black", ls=":", lw=1.0, label="e0 = %.4f" % e0)
ax2.set_xlabel("time t")
ax2.set_ylabel("energy per unit mass e")
ax2.set_ylim(e0 - 0.5, e0 + 0.5)  # zoom so the flat line is visually clear
ax2.legend(loc="upper right")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.1.1_s4.png")
print("saved figure to 5A.1.1_s4.png")
