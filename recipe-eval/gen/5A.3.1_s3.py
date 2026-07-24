import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Harmonic oscillator: d^2x/dt^2 = -k*x  <=>  dx/dt = v, dv/dt = -k*x
# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2  (constant for true motion)
# ----------------------------------------------------------------------

k  = 0.1        # spring stiffness
x0 = 1.0        # initial position
v0 = 2.0        # initial velocity
T  = 100.0      # final time

def accel(x):
    # acceleration f = dv/dt = -k*x
    return -k * x

def energy(x, v):
    # energy per unit mass
    return 0.5 * k * x**2 + 0.5 * v**2

# ----------------------------------------------------------------------
# Explicit (forward) Euler, for comparison
# ----------------------------------------------------------------------
def euler(dt):
    n = int(round(T / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1)
    v = np.empty(n + 1)
    x[0], v[0] = x0, v0
    for i in range(n):
        # update using values at the current step
        x[i + 1] = x[i] + dt * v[i]
        v[i + 1] = v[i] + dt * accel(x[i])
    return t, x, v

# ----------------------------------------------------------------------
# Leapfrog (velocity carried at half-integer steps)
#   startup half-step:  v_half = v0 + 0.5*dt*f0
#   then per step:       x_next        = x + dt*v_half
#                        v_next_half   = v_half + dt*f(x_next)
# We report x at integer steps, and v synced to integer steps by
# averaging the surrounding half-step velocities (only for the plotted
# "synced" energy); the raw staggered energy is also computed.
# ----------------------------------------------------------------------
def leapfrog(dt):
    n = int(round(T / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1)
    v_sync = np.empty(n + 1)     # velocity interpolated to integer steps
    e_stag = np.empty(n + 1)     # energy using staggered (half-step) velocity

    x[0] = x0
    v_prev_half = v0 + 0.5 * dt * accel(x0)   # startup half-step
    v_sync[0] = v0                            # known exactly at t=0
    e_stag[0] = energy(x0, v_prev_half)       # x at t=0, v at t=dt/2 (staggered)

    for i in range(n):
        # drift: advance position a full step using half-step velocity
        x[i + 1] = x[i] + dt * v_prev_half
        # kick: advance velocity a full step using new acceleration
        v_next_half = v_prev_half + dt * accel(x[i + 1])
        # velocity at integer step = average of the two surrounding half-steps
        v_sync[i + 1] = 0.5 * (v_prev_half + v_next_half)
        # staggered energy: x at integer step, v at the following half-step
        e_stag[i + 1] = energy(x[i + 1], v_next_half)
        v_prev_half = v_next_half

    return t, x, v_sync, e_stag

# ----------------------------------------------------------------------
# Run both methods at both step sizes
# ----------------------------------------------------------------------
e_exact = energy(x0, v0)
print(f"Exact conserved energy e0 = {e_exact:.6f}")

results = {}
for dt in (0.01, 0.1):
    tl, xl, vl, el_stag = leapfrog(dt)
    te, xe, ve = euler(dt)
    el_sync = energy(xl, vl)   # energy from synced (integer-step) x and v

    results[dt] = (tl, xl, el_sync, el_stag, te, xe, energy(xe, ve))

    print(f"--- dt = {dt} ---")
    print(f"leapfrog synced-energy  min = {el_sync.min():.6f}, max = {el_sync.max():.6f}, drift = {el_sync[-1]-el_sync[0]:+.3e}")
    print(f"leapfrog staggered-e    min = {el_stag.min():.6f}, max = {el_stag.max():.6f}, drift = {el_stag[-1]-el_stag[0]:+.3e}")
    print(f"euler energy            min = {energy(xe,ve).min():.6f}, max = {energy(xe,ve).max():.6f}, drift = {energy(xe,ve)[-1]-energy(xe,ve)[0]:+.3e}")

# Drift check: fit a straight line to leapfrog energy at the large step.
tl, xl, el_sync, el_stag, te, xe, ee = results[0.1]
slope = np.polyfit(tl, el_stag, 1)[0]
print(f"leapfrog (dt=0.1) staggered-energy linear drift slope = {slope:+.3e} per unit time (≈0 confirms no secular drift)")

# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

for col, dt in enumerate((0.01, 0.1)):
    tl, xl, el_sync, el_stag, te, xe, ee = results[dt]

    ax = axes[0, col]
    ax.plot(tl, xl, lw=1.0, label="leapfrog x(t)")
    ax.plot(te, xe, lw=0.8, alpha=0.7, label="euler x(t)")
    ax.set_title(f"Position x(t), dt = {dt}")
    ax.set_xlabel("t"); ax.set_ylabel("x"); ax.legend(); ax.grid(alpha=0.3)

    ax = axes[1, col]
    ax.plot(tl, el_stag, lw=1.0, label="leapfrog e(t) (staggered)")
    ax.plot(tl, el_sync, lw=1.0, label="leapfrog e(t) (synced)")
    ax.plot(te, ee, lw=0.8, alpha=0.7, label="euler e(t)")
    ax.axhline(e_exact, color="k", ls="--", lw=0.8, label="exact e0")
    ax.set_title(f"Energy e(t), dt = {dt}")
    ax.set_xlabel("t"); ax.set_ylabel("e"); ax.legend(); ax.grid(alpha=0.3)

fig.suptitle("Harmonic oscillator: leapfrog (symplectic) vs Euler energy behavior")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.3.1_s3.png", dpi=120)

# One-sentence explanation of why the check confirms the result:
print("Explanation: The leapfrog energy oscillates within a fixed bounded band with ~zero linear-drift slope even at dt=0.1, whereas Euler's energy grows without bound; this confirms leapfrog is symplectic/time-reversible and that the residual ripple is merely an artifact of x and v being stored at staggered times, not true energy loss.")
