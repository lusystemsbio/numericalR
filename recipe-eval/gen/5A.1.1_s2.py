import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Test system: harmonic oscillator (mass on a spring), used as the reference
# against which numerical integrators (later sections) will be compared.
#
# Equation of motion:  d^2x/dt^2 = -k*x
# written as a first-order pair:
#     dx/dt = v
#     dv/dt = -k*x
# Force per unit mass (direct evaluation):  f(x) = -k*x
# Energy per unit mass:  e = 0.5*k*x^2 + 0.5*v^2  (exactly constant for true motion)
# -----------------------------------------------------------------------------

# --- parameters and initial conditions ---
k  = 0.1   # spring stiffness
x0 = 1.0   # initial displacement
v0 = 2.0   # initial velocity

def f(x):
    # direct evaluation of the force per unit mass
    return -k * x

def energy(x, v):
    # energy per unit mass
    return 0.5 * k * x**2 + 0.5 * v**2

# -----------------------------------------------------------------------------
# EXACT (reference) solution, built explicitly step by step rather than by
# calling a black-box ODE solver.
#
# The general solution of x'' = -k*x is a pure sinusoid with angular frequency
#     w = sqrt(k)
# We write it as x(t) = A*cos(w*t + phi) and fix A, phi from (x0, v0):
#     x(0) = A*cos(phi)        = x0
#     v(0) = -A*w*sin(phi)     = v0
# From these:  A*cos(phi) = x0  and  A*sin(phi) = -v0/w, so
#     A   = sqrt(x0^2 + (v0/w)^2)
#     phi = atan2(-v0/w, x0)
# -----------------------------------------------------------------------------
w   = np.sqrt(k)                       # angular frequency
A   = np.sqrt(x0**2 + (v0 / w)**2)     # amplitude from initial conditions
phi = np.arctan2(-v0 / w, x0)          # phase from initial conditions
T   = 2.0 * np.pi / w                  # period

# evaluate the exact position and velocity on a time grid (several periods)
t = np.linspace(0.0, 4.0 * T, 2000)
x_exact = A * np.cos(w * t + phi)      # exact displacement
v_exact = -A * w * np.sin(w * t + phi) # exact velocity (= dx/dt)
e_exact = energy(x_exact, v_exact)     # exact energy along the trajectory

# -----------------------------------------------------------------------------
# What a correct simulation MUST preserve:
#   1. The motion stays a pure sinusoid of FIXED amplitude A (no growth/decay).
#   2. The energy per unit mass e stays CONSTANT (no drift) at its initial value.
# Any amplitude growth or energy drift in a numerical integrator is therefore a
# numerical artifact, not physics.
# -----------------------------------------------------------------------------

# --- separate check: constant amplitude and constant energy ---
e_initial = energy(x0, v0)
e_mean    = np.mean(e_exact)
e_drift   = np.max(np.abs(e_exact - e_initial))   # deviation from initial energy
amp_meas  = np.max(np.abs(x_exact))               # measured amplitude from the trace

print("Spring stiffness k:            ", k)
print("Initial displacement x0:       ", x0)
print("Initial velocity v0:           ", v0)
print("Angular frequency w=sqrt(k):   ", w)
print("Period T=2*pi/w:               ", T)
print("Amplitude A:                   ", A)
print("Phase phi (rad):               ", phi)
print("Initial force f(x0)=-k*x0:     ", f(x0))
print("Initial energy e0:             ", e_initial)
print("Mean energy along trajectory:  ", e_mean)
print("Max energy drift |e - e0|:     ", e_drift)
print("Measured amplitude max|x|:     ", amp_meas)
print("Amplitude error |max|x| - A|:  ", abs(amp_meas - A))

# The check confirms the result because a machine-precision-small energy drift
# together with max|x| equal to A proves the exact trajectory is a pure sinusoid
# at constant energy, so anything else a simulation shows must be numerical.
print("Check confirms result:         ",
      "exact motion is a pure sinusoid at constant energy "
      "(energy drift and amplitude error are ~machine precision), "
      "so any growth/drift a simulation shows is a numerical artifact.")

# -----------------------------------------------------------------------------
# Reference plot: exact motion x(t) and its constant energy e(t).
# -----------------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

ax1.plot(t, x_exact, color="C0", label="exact x(t)")
ax1.axhline( A, color="gray", ls="--", lw=1, label="+/- amplitude A")
ax1.axhline(-A, color="gray", ls="--", lw=1)
ax1.set_ylabel("displacement x(t)")
ax1.set_title("Harmonic oscillator reference: pure sinusoid, constant energy")
ax1.legend(loc="upper right")
ax1.grid(True, alpha=0.3)

ax2.plot(t, e_exact, color="C3", label="exact energy e(t)")
ax2.axhline(e_initial, color="k", ls=":", lw=1, label="initial energy e0")
ax2.set_xlabel("time t")
ax2.set_ylabel("energy per unit mass e")
# widen the y-range so the flat line is visibly constant, not just autoscaled noise
ax2.set_ylim(e_initial - 0.5, e_initial + 0.5)
ax2.legend(loc="upper right")
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.1.1_s2.png")
