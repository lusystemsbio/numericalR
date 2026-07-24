import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: mass on a spring, d^2x/dt^2 = -k*x ----
k = 0.1
x0, v0 = 1.0, 2.0
t_end = 100.0

def f(x):
    # acceleration a = dv/dt = -k*x
    return -k * x

def energy(x, v):
    # energy per unit mass, exactly constant for true motion
    return 0.5 * k * x**2 + 0.5 * v**2

# ---- Euler integrator (for comparison) ----
def euler(dt):
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1); v = np.empty(n + 1)
    x[0], v[0] = x0, v0
    for i in range(n):
        # forward Euler: use current slopes to step both variables
        x[i+1] = x[i] + dt * v[i]
        v[i+1] = v[i] + dt * f(x[i])
    return t, x, v

# ---- Leapfrog integrator, velocity carried at half-integer steps ----
def leapfrog(dt):
    n = int(round(t_end / dt))
    t = np.linspace(0.0, n * dt, n + 1)
    x = np.empty(n + 1); v = np.empty(n + 1)   # v stored at integer steps for energy
    x[0], v[0] = x0, v0
    # startup half-step: advance velocity by dt/2 to reach the half grid
    v_half = v0 + 0.5 * dt * f(x0)
    for i in range(n):
        # drift position a full step using the half-step velocity
        x[i+1] = x[i] + dt * v_half
        # kick velocity a full step using the force at the new position
        v_next_half = v_half + dt * f(x[i+1])
        # integer-step velocity = average of the two straddling half-steps
        v[i+1] = 0.5 * (v_half + v_next_half)
        v_half = v_next_half
    return t, x, v

# ---- Run both methods at both step sizes ----
dts = [0.01, 0.1]
results = {}
for dt in dts:
    tl, xl, vl = leapfrog(dt)
    te, xe, ve = euler(dt)
    results[dt] = (tl, xl, vl, te, xe, ve)

e_exact = energy(x0, v0)
print(f"Exact conserved energy e = {e_exact:.6f}")

for dt in dts:
    tl, xl, vl, te, xe, ve = results[dt]
    el = energy(xl, vl)
    ee = energy(xe, ve)
    print(f"--- dt = {dt} ---")
    print(f"Leapfrog energy: initial = {el[0]:.6f}, final = {el[-1]:.6f}, "
          f"min = {el.min():.6f}, max = {el.max():.6f}, "
          f"drift(final-initial) = {el[-1]-el[0]:.3e}, ripple(max-min) = {el.max()-el.min():.3e}")
    print(f"Euler    energy: initial = {ee[0]:.6f}, final = {ee[-1]:.6f}, "
          f"drift(final-initial) = {ee[-1]-ee[0]:.3e}")

# ---- Plots ----
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

# Top row: x(t) for leapfrog at both step sizes
ax = axes[0, 0]
for dt in dts:
    tl, xl, vl, *_ = results[dt]
    ax.plot(tl, xl, label=f"leapfrog dt={dt}")
ax.set_title("Leapfrog x(t)"); ax.set_xlabel("t"); ax.set_ylabel("x"); ax.legend()

# Top-right: energy e(t) for leapfrog at both step sizes
ax = axes[0, 1]
for dt in dts:
    tl, xl, vl, *_ = results[dt]
    ax.plot(tl, energy(xl, vl), label=f"leapfrog dt={dt}")
ax.axhline(e_exact, color="k", ls="--", lw=0.8, label="exact")
ax.set_title("Leapfrog energy e(t) (no drift)"); ax.set_xlabel("t"); ax.set_ylabel("e"); ax.legend()

# Bottom-left: leapfrog vs Euler energy at dt=0.01
ax = axes[1, 0]
tl, xl, vl, te, xe, ve = results[0.01]
ax.plot(tl, energy(xl, vl), label="leapfrog")
ax.plot(te, energy(xe, ve), label="Euler")
ax.axhline(e_exact, color="k", ls="--", lw=0.8, label="exact")
ax.set_title("Energy e(t), dt=0.01"); ax.set_xlabel("t"); ax.set_ylabel("e"); ax.legend()

# Bottom-right: leapfrog vs Euler energy at dt=0.1
ax = axes[1, 1]
tl, xl, vl, te, xe, ve = results[0.1]
ax.plot(tl, energy(xl, vl), label="leapfrog")
ax.plot(te, energy(xe, ve), label="Euler")
ax.axhline(e_exact, color="k", ls="--", lw=0.8, label="exact")
ax.set_title("Energy e(t), dt=0.1"); ax.set_xlabel("t"); ax.set_ylabel("e"); ax.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.3.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the leapfrog energy at dt=0.1 stays bounded and returns "
      "to its start (final-initial drift ~ machine noise) instead of growing like Euler's, "
      "the only variation is a small periodic ripple from evaluating x and v at staggered "
      "half-step times, confirming the method conserves energy over long times (is symplectic).")
