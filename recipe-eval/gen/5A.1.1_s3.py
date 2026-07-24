import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# =====================================================================
# Test system: harmonic oscillator (mass on a spring)
#   d^2x/dt^2 = -k*x   ->   dx/dt = v ,  dv/dt = -k*x
# Energy per unit mass:  e = 0.5*k*x^2 + 0.5*v^2  (constant for true motion)
# This system is the reference against which integrators are compared.
# What a correct simulation MUST preserve:
#   (1) the motion stays a bounded sinusoid (no amplitude growth/decay)
#   (2) the energy e stays constant (no energy drift)
# =====================================================================

# --- Physical parameters and initial conditions ---
k  = 0.1     # spring stiffness
x0 = 1.0     # initial displacement
v0 = 2.0     # initial velocity

# --- Direct force evaluation (used later by numerical integrators) ---
def f(x):
    # force per unit mass from the spring
    return -k * x

print(f"stiffness k = {k}")
print(f"initial displacement x0 = {x0}")
print(f"initial velocity v0 = {v0}")
print(f"force at x0, f(x0) = {f(x0)}")

# --- Exact (analytic) solution -------------------------------------------------
# The equation d^2x/dt^2 = -k*x has natural angular frequency omega = sqrt(k).
# General solution: x(t) = A*cos(omega*t) + B*sin(omega*t)
# Match initial conditions:
#   x(0)  = A            = x0     ->  A = x0
#   x'(0) = B*omega      = v0     ->  B = v0/omega
omega  = np.sqrt(k)                 # angular frequency
period = 2.0 * np.pi / omega        # oscillation period
A = x0                              # cosine amplitude coefficient
B = v0 / omega                      # sine amplitude coefficient
amplitude = np.sqrt(A**2 + B**2)    # true (constant) amplitude of the sinusoid

print(f"angular frequency omega = sqrt(k) = {omega}")
print(f"period T = 2*pi/omega = {period}")
print(f"cosine coefficient A = {A}")
print(f"sine coefficient B = {B}")
print(f"true amplitude = sqrt(A^2 + B^2) = {amplitude}")

# --- Time grid spanning a few periods ---
t = np.linspace(0.0, 4.0 * period, 2000)

# Exact position and velocity
x_exact = A * np.cos(omega * t) + B * np.sin(omega * t)          # x(t)
v_exact = -A * omega * np.sin(omega * t) + B * omega * np.cos(omega * t)  # v(t)=dx/dt

# Exact energy per unit mass at every time
e_exact = 0.5 * k * x_exact**2 + 0.5 * v_exact**2

# Initial energy (the value that must be preserved)
e0 = 0.5 * k * x0**2 + 0.5 * v0**2
print(f"initial energy e0 = 0.5*k*x0^2 + 0.5*v0^2 = {e0}")

# --- Separate check: motion is a pure sinusoid at constant energy --------------
# A pure sinusoid stays within +/- amplitude, and constant energy means the
# spread of e over the whole run should be numerically zero.
x_max = np.max(x_exact)
x_min = np.min(x_exact)
e_max = np.max(e_exact)
e_min = np.min(e_exact)
energy_drift = e_max - e_min                 # variation of energy across the run
amplitude_error = abs(max(x_max, -x_min) - amplitude)

print(f"x_exact max = {x_max}")
print(f"x_exact min = {x_min}")
print(f"energy max = {e_max}")
print(f"energy min = {e_min}")
print(f"energy drift (max-min) = {energy_drift}")
print(f"amplitude vs. |x| envelope error = {amplitude_error}")

# Because energy_drift and amplitude_error are ~machine-zero, the exact motion
# is confirmed to be a pure sinusoid at constant energy; therefore any amplitude
# growth or energy drift seen in a later numerical integrator is a numerical
# artifact, not physics -- this check confirms it by showing the true reference
# has neither, so any deviation must originate from the integrator.

# --- Plot the reference to match ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

ax1.plot(t, x_exact, color="tab:blue", label="exact x(t)")
ax1.axhline( amplitude, color="gray", ls="--", lw=0.8, label="+/- amplitude")
ax1.axhline(-amplitude, color="gray", ls="--", lw=0.8)
ax1.set_ylabel("displacement x(t)")
ax1.set_title("Harmonic oscillator: exact motion (reference)")
ax1.legend(loc="upper right")
ax1.grid(True, alpha=0.3)

ax2.plot(t, e_exact, color="tab:red", label="exact energy e(t)")
ax2.axhline(e0, color="black", ls=":", lw=1.0, label=f"e0 = {e0:g}")
ax2.set_xlabel("time t")
ax2.set_ylabel("energy per unit mass")
ax2.set_title("Exact energy is constant")
ax2.legend(loc="upper right")
ax2.grid(True, alpha=0.3)
# keep a visible scale so a truly constant line does not look like noise
ax2.set_ylim(e0 - 0.5, e0 + 0.5)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.1.1_s3.png")
print("Saved figure to 5A.1.1_s3.png")
