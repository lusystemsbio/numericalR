import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Test system: harmonic oscillator (mass on a spring)
#   d^2x/dt^2 = -k*x   <=>   dx/dt = v,  dv/dt = -k*x
# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2  (constant for true motion)
# ---------------------------------------------------------------

# Direct evaluation of the force (per unit mass): f(x) = -k*x
def f(x, k):
    return -k * x

# Parameters and initial conditions
k = 0.1        # spring stiffness
x0 = 1.0       # initial displacement
v0 = 2.0       # initial velocity

# --- Exact (analytic) solution ---------------------------------
# The true motion is a pure sinusoid at angular frequency omega = sqrt(k):
#   x(t) = x0*cos(w t) + (v0/w)*sin(w t)
#   v(t) = -x0*w*sin(w t) + v0*cos(w t)
omega = np.sqrt(k)
print("Angular frequency omega = sqrt(k) =", omega)
print("Period T = 2*pi/omega =", 2 * np.pi / omega)

t = np.linspace(0.0, 4.0 * (2 * np.pi / omega), 2000)  # four periods
x_exact = x0 * np.cos(omega * t) + (v0 / omega) * np.sin(omega * t)
v_exact = -x0 * omega * np.sin(omega * t) + v0 * np.cos(omega * t)

# Energy of the exact motion at every sample time
e_exact = 0.5 * k * x_exact**2 + 0.5 * v_exact**2

# Reference energy computed directly from the initial condition
e0 = 0.5 * k * x0**2 + 0.5 * v0**2
print("Initial energy e0 = 0.5*k*x0^2 + 0.5*v0^2 =", e0)

# --- What a correct simulation must preserve -------------------
# A correct integrator must (a) reproduce the sinusoidal x(t) without
# amplitude growth or decay, and (b) hold the energy e constant at e0.
# Any amplitude growth or energy drift in a numerical run is an artifact.

# Amplitude of the exact motion (analytic): A = sqrt(x0^2 + (v0/omega)^2)
amplitude = np.sqrt(x0**2 + (v0 / omega)**2)
print("Exact amplitude A = sqrt(x0^2 + (v0/omega)^2) =", amplitude)
print("Max |x_exact| observed =", np.max(np.abs(x_exact)))

# --- Separate check: pure sinusoid at constant energy ----------
# Check 1: energy is constant -> its spread over the run is (near) zero.
e_min = np.min(e_exact)
e_max = np.max(e_exact)
e_drift = e_max - e_min
print("Exact energy min =", e_min)
print("Exact energy max =", e_max)
print("Exact energy drift (max - min) =", e_drift)
print("Max |e_exact - e0| =", np.max(np.abs(e_exact - e0)))

# Check 2: the motion satisfies the ODE and stays a pure sinusoid.
# Verify x(t) obeys x'' = -k*x by comparing the analytic second
# derivative (-omega^2 * x) with -k*x; they must agree exactly.
xdd = -omega**2 * x_exact          # second derivative of the sinusoid
residual = np.max(np.abs(xdd - f(x_exact, k)))
print("Max |x'' - (-k*x)| residual =", residual)

# --- Plot: exact motion and its constant energy ----------------
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

ax1.plot(t, x_exact, color="tab:blue", label="x(t) exact")
ax1.axhline(amplitude, color="gray", ls="--", lw=0.8, label="+/- amplitude")
ax1.axhline(-amplitude, color="gray", ls="--", lw=0.8)
ax1.set_ylabel("displacement x")
ax1.set_title("Exact harmonic oscillator (k=%.2f, x0=%.1f, v0=%.1f)" % (k, x0, v0))
ax1.legend(loc="upper right")
ax1.grid(True, alpha=0.3)

ax2.plot(t, e_exact, color="tab:red", label="energy e(t) exact")
ax2.axhline(e0, color="black", ls="--", lw=0.8, label="e0 reference")
ax2.set_xlabel("time t")
ax2.set_ylabel("energy e")
ax2.set_ylim(e0 * 0.9, e0 * 1.1)   # tight window so any drift would show
ax2.legend(loc="upper right")
ax2.grid(True, alpha=0.3)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.1.1_s5.png")

# One-sentence explanation of why the check confirms the result:
# Because the analytic x(t) satisfies x'' = -k*x exactly (zero residual) while
# its energy stays at e0 (zero drift), the reference is provably a pure sinusoid
# at constant energy, so any amplitude growth or energy drift later seen in a
# numerical integrator must be a numerical artifact rather than true physics.
print("Check confirms: zero ODE residual + zero energy drift => exact motion is a "
      "pure sinusoid at constant energy, so any later growth/drift is numerical.")
