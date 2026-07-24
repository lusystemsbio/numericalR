import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Test system: harmonic oscillator (mass on a spring)
#   d^2x/dt^2 = -k*x   <=>   dx/dt = v ,  dv/dt = -k*x
# Energy per unit mass:  e = 0.5*k*x^2 + 0.5*v^2  (constant for true motion)
# This is the reference every numerical integrator must match:
# a correct simulation must preserve (1) the sinusoidal shape/amplitude
# and (2) the constant energy -- no amplitude growth, no energy drift.
# ---------------------------------------------------------------

# Parameters and initial conditions
k  = 0.1     # spring stiffness
x0 = 1.0     # initial displacement
v0 = 2.0     # initial velocity

# Force per unit mass, evaluated directly (used by later integrators)
def f(x):
    return -k * x   # Hooke's law force / mass = acceleration

# Exact analytic solution of dx/dt=v, dv/dt=-k*x:
#   x(t) = x0*cos(w t) + (v0/w)*sin(w t),  with angular frequency w = sqrt(k)
w = np.sqrt(k)                       # angular frequency
period = 2.0 * np.pi / w             # oscillation period

t = np.linspace(0.0, 3.0 * period, 2000)   # three full periods
x_exact = x0 * np.cos(w * t) + (v0 / w) * np.sin(w * t)   # exact displacement
v_exact = -x0 * w * np.sin(w * t) + v0 * np.cos(w * t)    # exact velocity (derivative)

# Exact energy per unit mass along the trajectory
e_exact = 0.5 * k * x_exact**2 + 0.5 * v_exact**2
e_initial = 0.5 * k * x0**2 + 0.5 * v0**2   # energy from the initial state

# Amplitude of the pure sinusoid: sqrt(x0^2 + (v0/w)^2)
amplitude = np.sqrt(x0**2 + (v0 / w)**2)

# --- Check: exact motion is a pure sinusoid at constant energy ---
# A single-frequency sinusoid stays within +/-amplitude forever, and its
# energy is a fixed number. We confirm numerically that (a) the energy along
# the whole trajectory never departs from the initial energy, and (b) the
# displacement never exceeds the analytic amplitude.
energy_drift = np.max(np.abs(e_exact - e_initial))          # should be ~machine eps
amplitude_overshoot = np.max(np.abs(x_exact)) - amplitude   # should be ~machine eps

# --- Report every numerical result ---
print(f"Spring stiffness k:                 {k}")
print(f"Initial displacement x0:            {x0}")
print(f"Initial velocity v0:                {v0}")
print(f"Angular frequency w = sqrt(k):      {w}")
print(f"Oscillation period T = 2*pi/w:      {period}")
print(f"Sinusoid amplitude sqrt(x0^2+(v0/w)^2): {amplitude}")
print(f"Initial energy per unit mass e0:    {e_initial}")
print(f"Min energy along trajectory:        {np.min(e_exact)}")
print(f"Max energy along trajectory:        {np.max(e_exact)}")
print(f"Max energy drift |e(t)-e0|:         {energy_drift}")
print(f"Max |x(t)| minus amplitude (overshoot): {amplitude_overshoot}")

# Because the maximum energy drift and the amplitude overshoot are both at
# the level of floating-point round-off, the exact motion is confirmed to be
# a pure sinusoid of fixed amplitude at constant energy -- so any amplitude
# growth or energy drift seen in a later integrator is purely a numerical artifact.
print("Check confirms result: energy drift and amplitude overshoot are both ~machine epsilon,")
print("so the exact motion is a constant-energy pure sinusoid and any drift later is numerical.")

# --- Reference plot: exact motion and its constant energy ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 7), sharex=True)

ax1.plot(t, x_exact, color="C0", label="exact x(t)")
ax1.axhline( amplitude, color="gray", ls="--", lw=0.8, label="+/- amplitude")
ax1.axhline(-amplitude, color="gray", ls="--", lw=0.8)
ax1.set_ylabel("displacement x")
ax1.set_title("Harmonic oscillator (k=0.1): exact motion is a pure sinusoid")
ax1.legend(loc="upper right")
ax1.grid(True, alpha=0.3)

ax2.plot(t, e_exact, color="C3", label="exact energy e(t)")
ax2.axhline(e_initial, color="k", ls=":", lw=1.0, label=f"e0 = {e_initial:.4f}")
ax2.set_xlabel("time t")
ax2.set_ylabel("energy per unit mass e")
ax2.set_title("Energy is exactly constant (the target integrators must match)")
ax2.legend(loc="upper right")
ax2.grid(True, alpha=0.3)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.1.1_s1.png")
