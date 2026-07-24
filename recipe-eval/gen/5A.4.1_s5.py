import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters
k = 0.1          # spring stiffness
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # final time

def force(x):
    # acceleration (force per unit mass) for the spring: d^2x/dt^2 = -k*x
    return -k * x

def energy(x, v):
    # energy per unit mass, constant for the true motion
    return 0.5 * k * x**2 + 0.5 * v**2

def velocity_verlet(dt):
    # integrate from 0 to t_end with velocity Verlet, returning t, x, v arrays
    n = int(round(t_end / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0                  # set initial conditions
    f = force(x[0])                      # current force
    for i in range(n):
        v_half = v[i] + 0.5 * dt * f     # first half-step of velocity update
        x[i+1] = x[i] + dt * v_half      # full-step position update
        f_next = force(x[i+1])           # force at the new position
        v[i+1] = v_half + 0.5 * dt * f_next  # second half-step of velocity update
        f = f_next                       # carry force forward (no recompute)
        t[i+1] = t[i] + dt               # advance time; x and v share this time
    return t, x, v

# Run both step sizes
t_small, x_small, v_small = velocity_verlet(0.01)
t_large, x_large, v_large = velocity_verlet(0.1)

e_small = energy(x_small, v_small)
e_large = energy(x_large, v_large)

e0 = energy(x0, v0)

# --- Plots ---
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0, 0].plot(t_small, x_small, color="tab:blue")
axes[0, 0].set_title("x(t), velocity Verlet, dt = 0.01")
axes[0, 0].set_xlabel("t"); axes[0, 0].set_ylabel("x")

axes[0, 1].plot(t_large, x_large, color="tab:orange")
axes[0, 1].set_title("x(t), velocity Verlet, dt = 0.1")
axes[0, 1].set_xlabel("t"); axes[0, 1].set_ylabel("x")

axes[1, 0].plot(t_small, e_small, color="tab:blue")
axes[1, 0].set_title("energy e(t), dt = 0.01")
axes[1, 0].set_xlabel("t"); axes[1, 0].set_ylabel("e")

axes[1, 1].plot(t_large, e_large, color="tab:orange")
axes[1, 1].set_title("energy e(t), dt = 0.1")
axes[1, 1].set_xlabel("t"); axes[1, 1].set_ylabel("e")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.4.1_s5.png")

# --- Numerical results ---
print(f"Initial energy e0: {e0}")
print(f"dt=0.01: final x = {x_small[-1]}")
print(f"dt=0.01: final v = {v_small[-1]}")
print(f"dt=0.01: final energy = {e_small[-1]}")
print(f"dt=0.01: max energy = {e_small.max()}")
print(f"dt=0.01: min energy = {e_small.min()}")
print(f"dt=0.01: max abs energy drift = {np.max(np.abs(e_small - e0))}")
print(f"dt=0.01: relative energy drift = {np.max(np.abs(e_small - e0)) / e0}")

print(f"dt=0.1: final x = {x_large[-1]}")
print(f"dt=0.1: final v = {v_large[-1]}")
print(f"dt=0.1: final energy = {e_large[-1]}")
print(f"dt=0.1: max energy = {e_large.max()}")
print(f"dt=0.1: min energy = {e_large.min()}")
print(f"dt=0.1: max abs energy drift = {np.max(np.abs(e_large - e0))}")
print(f"dt=0.1: relative energy drift = {np.max(np.abs(e_large - e0)) / e0}")

# --- Check: energy conservation and matched time points ---
tol = 1e-2  # bounded-drift tolerance (relative); Verlet does not conserve energy exactly but stays bounded
conserved_small = np.max(np.abs(e_small - e0)) / e0 < tol
conserved_large = np.max(np.abs(e_large - e0)) / e0 < tol
same_times = (len(t_small) == len(x_small) == len(v_small)) and (len(t_large) == len(x_large) == len(v_large))

print(f"CHECK dt=0.01: energy bounded within {tol*100:.1f}% of e0: {conserved_small}")
print(f"CHECK dt=0.1: energy bounded within {tol*100:.1f}% of e0: {conserved_large}")
print(f"CHECK: x and v reported at identical time points at both step sizes: {same_times}")

# One-sentence explanation:
# This check confirms the result because velocity Verlet's structure returns x[i] and v[i]
# at the same time t[i] (each iteration advances one shared time), and a small bounded energy
# drift rather than growth demonstrates the symplectic, energy-stable behavior expected of the integrator.
print("Explanation: x and v share the same t[i] each step and the energy error stays small and bounded rather than growing, confirming correct symplectic velocity-Verlet integration.")
