import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters: mass on a spring, d^2x/dt^2 = -k*x
k = 0.1          # spring stiffness
x0 = 1.0         # initial position
v0 = 2.0         # initial velocity
t_end = 100.0    # final time

# Energy per unit mass e = 0.5*k*x^2 + 0.5*v^2 (exactly constant for the true motion)
def energy(x, v):
    return 0.5 * k * x**2 + 0.5 * v**2

# Explicit (forward) Euler integrator for Newton's equations.
# We update velocity first using the force f = -k*x at the current position,
# then update position using the OLD velocity, matching the requested scheme:
#   v_next = v + dt*f
#   x_next = x + dt*v
def euler_integrate(dt):
    n = int(round(t_end / dt))          # number of steps
    t = np.zeros(n + 1)
    x = np.zeros(n + 1)
    v = np.zeros(n + 1)
    x[0], v[0] = x0, v0                  # set initial conditions
    for i in range(n):
        f = -k * x[i]                   # force per unit mass at current position
        v[i + 1] = v[i] + dt * f        # 1) update velocity from the force
        x[i + 1] = x[i] + dt * v[i]     # 2) update position from the OLD velocity
        t[i + 1] = t[i] + dt            # advance time
    return t, x, v

# Run both step sizes
t_fine, x_fine, v_fine = euler_integrate(0.01)
t_coarse, x_coarse, v_coarse = euler_integrate(0.1)

e_fine = energy(x_fine, v_fine)
e_coarse = energy(x_coarse, v_coarse)

e_initial = energy(x0, v0)

# Report the energy check numerically
print(f"Initial energy e(0)                 = {e_initial:.6f}")
print(f"dt=0.01: final energy e(t_end)      = {e_fine[-1]:.6f}")
print(f"dt=0.01: energy change (final-init) = {e_fine[-1] - e_initial:.6f}")
print(f"dt=0.01: max |x| over run           = {np.max(np.abs(x_fine)):.6f}")
print(f"dt=0.10: final energy e(t_end)      = {e_coarse[-1]:.6f}")
print(f"dt=0.10: energy change (final-init) = {e_coarse[-1] - e_initial:.6f}")
print(f"dt=0.10: max |x| over run           = {np.max(np.abs(x_coarse)):.6f}")

# Growth factors make the non-conservation explicit
print(f"dt=0.01: energy growth factor e_final/e_init = {e_fine[-1] / e_initial:.6f}")
print(f"dt=0.10: energy growth factor e_final/e_init = {e_coarse[-1] / e_initial:.6f}")

# Plots: x(t) and e(t) for both step sizes
fig, axes = plt.subplots(2, 2, figsize=(12, 8))

axes[0, 0].plot(t_fine, x_fine, lw=0.8)
axes[0, 0].set_title("x(t), dt = 0.01")
axes[0, 0].set_xlabel("t")
axes[0, 0].set_ylabel("x")

axes[0, 1].plot(t_coarse, x_coarse, lw=0.8, color="tab:red")
axes[0, 1].set_title("x(t), dt = 0.1 (amplitude grows)")
axes[0, 1].set_xlabel("t")
axes[0, 1].set_ylabel("x")

axes[1, 0].plot(t_fine, e_fine, lw=0.8)
axes[1, 0].axhline(e_initial, color="k", ls="--", lw=0.8, label="true energy")
axes[1, 0].set_title("e(t), dt = 0.01 (slow upward drift)")
axes[1, 0].set_xlabel("t")
axes[1, 0].set_ylabel("e")
axes[1, 0].legend()

axes[1, 1].plot(t_coarse, e_coarse, lw=0.8, color="tab:red")
axes[1, 1].axhline(e_initial, color="k", ls="--", lw=0.8, label="true energy")
axes[1, 1].set_title("e(t), dt = 0.1 (energy grows steadily)")
axes[1, 1].set_xlabel("t")
axes[1, 1].set_ylabel("e")
axes[1, 1].legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5A.2.1_s5.png")

# One-sentence explanation of why the check confirms the result:
# Because the exact motion keeps e constant, any steady growth of the computed
# energy e(t) above its initial value can only come from the integrator itself,
# so seeing e drift/grow (fast at dt=0.1, slowly at dt=0.01) proves that explicit
# Euler injects spurious energy and does not conserve it.
print("Explanation: since true e is exactly constant, the observed upward drift/growth "
      "of e(t) (rapid at dt=0.1, slow at dt=0.01) can only be numerical error from the "
      "integrator, confirming that explicit Euler does not conserve energy.")
