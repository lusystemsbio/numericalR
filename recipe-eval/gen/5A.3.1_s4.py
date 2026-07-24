import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------
# Harmonic oscillator: dx/dt = v, dv/dt = -k*x
# Energy per unit mass: e = 0.5*k*x^2 + 0.5*v^2 (constant for true motion)
# ------------------------------------------------------------------
k = 0.1
x0, v0 = 1.0, 2.0
t_end = 100.0

def force(x):
    # acceleration f = -k*x
    return -k * x

def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

# ------------------------------------------------------------------
# Explicit (forward) Euler integrator -- for comparison
# ------------------------------------------------------------------
def euler(dt):
    n = int(round(t_end / dt))
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0
    for i in range(n):
        # update both using values at the current (integer) time
        x[i + 1] = x[i] + dt * v[i]          # x_{n+1} = x_n + dt*v_n
        v[i + 1] = v[i] + dt * force(x[i])   # v_{n+1} = v_n + dt*f_n
        t[i + 1] = t[i] + dt
    return t, x, v, energy(x, v)

# ------------------------------------------------------------------
# Leapfrog integrator -- velocity carried at half-integer steps
#   startup half-step: v_half = v0 + 0.5*dt*f0
#   then repeatedly:   x_next        = x + dt*v_half
#                      v_next_half   = v_half + dt*f_next
# ------------------------------------------------------------------
def leapfrog(dt):
    n = int(round(t_end / dt))
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)        # velocity stored at half-steps (staggered in time)
    x[0], v[0] = x0, v0
    # startup half-step to prime the staggered velocity
    v_half = v0 + 0.5 * dt * force(x0)
    for i in range(n):
        x[i + 1] = x[i] + dt * v_half            # advance position a full step
        v_next_half = v_half + dt * force(x[i + 1])  # advance half-step velocity
        # store a whole-step velocity estimate (average of the two half-steps)
        v[i + 1] = 0.5 * (v_half + v_next_half)
        v_half = v_next_half
        t[i + 1] = t[i] + dt
    return t, x, v, energy(x, v)

# ------------------------------------------------------------------
# Run for the two step sizes
# ------------------------------------------------------------------
tl_s, xl_s, vl_s, el_s = leapfrog(0.01)   # small step
tl_l, xl_l, vl_l, el_l = leapfrog(0.1)    # large step
te_s, xe_s, ve_s, ee_s = euler(0.01)
te_l, xe_l, ve_l, ee_l = euler(0.1)

e_exact = energy(x0, v0)

# ------------------------------------------------------------------
# Numerical diagnostics
# ------------------------------------------------------------------
print(f"Exact conserved energy e0                       : {e_exact:.6f}")
print(f"Leapfrog dt=0.01  mean energy                   : {np.mean(el_s):.6f}")
print(f"Leapfrog dt=0.01  energy drift (last - first)   : {el_s[-1] - el_s[0]:.6e}")
print(f"Leapfrog dt=0.01  max |e - e0|                  : {np.max(np.abs(el_s - e_exact)):.6e}")
print(f"Leapfrog dt=0.1   mean energy                   : {np.mean(el_l):.6f}")
print(f"Leapfrog dt=0.1   energy drift (last - first)   : {el_l[-1] - el_l[0]:.6e}")
print(f"Leapfrog dt=0.1   max |e - e0| (ripple only)    : {np.max(np.abs(el_l - e_exact)):.6e}")
print(f"Euler    dt=0.01  energy drift (last - first)   : {ee_s[-1] - ee_s[0]:.6e}")
print(f"Euler    dt=0.01  final energy                  : {ee_s[-1]:.6f}")
print(f"Euler    dt=0.1   energy drift (last - first)   : {ee_l[-1] - ee_l[0]:.6e}")
print(f"Euler    dt=0.1   final energy                  : {ee_l[-1]:.6f}")

# Confirm leapfrog does not drift secularly even at the large step:
# fit a straight line to e(t); the slope measures long-term drift per unit time.
slope_lf_large = np.polyfit(tl_l, el_l, 1)[0]
slope_eu_large = np.polyfit(te_l, ee_l, 1)[0]
print(f"Leapfrog dt=0.1   energy trend (slope de/dt)    : {slope_lf_large:.6e}")
print(f"Euler    dt=0.1   energy trend (slope de/dt)    : {slope_eu_large:.6e}")
# Explanation: the leapfrog slope is ~0 (bounded oscillation) while Euler's is
# strongly positive, so a near-zero slope with only a small bounded ripple
# confirms leapfrog conserves energy over time rather than drifting.

# ------------------------------------------------------------------
# Plots
# ------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

axes[0, 0].plot(tl_s, xl_s, lw=0.8)
axes[0, 0].set_title("Leapfrog x(t), dt=0.01")
axes[0, 0].set_xlabel("t"); axes[0, 0].set_ylabel("x")

axes[0, 1].plot(tl_l, xl_l, lw=0.8, color="tab:orange")
axes[0, 1].set_title("Leapfrog x(t), dt=0.1")
axes[0, 1].set_xlabel("t"); axes[0, 1].set_ylabel("x")

axes[1, 0].plot(tl_s, el_s, lw=0.8, label="leapfrog dt=0.01")
axes[1, 0].plot(tl_l, el_l, lw=0.8, label="leapfrog dt=0.1")
axes[1, 0].axhline(e_exact, color="k", ls="--", lw=0.8, label="exact e0")
axes[1, 0].set_title("Leapfrog energy e(t) -- no drift")
axes[1, 0].set_xlabel("t"); axes[1, 0].set_ylabel("e")
axes[1, 0].legend(fontsize=8)

axes[1, 1].plot(te_s, ee_s, lw=0.8, label="Euler dt=0.01")
axes[1, 1].plot(te_l, ee_l, lw=0.8, label="Euler dt=0.1")
axes[1, 1].plot(tl_l, el_l, lw=0.8, label="leapfrog dt=0.1")
axes[1, 1].axhline(e_exact, color="k", ls="--", lw=0.8, label="exact e0")
axes[1, 1].set_title("Euler energy grows vs. leapfrog stays bounded")
axes[1, 1].set_xlabel("t"); axes[1, 1].set_ylabel("e")
axes[1, 1].legend(fontsize=8)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.3.1_s4.png", dpi=120)
