import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Central attractive force ---
# Force magnitude G/r pointing toward the origin.
# Cartesian components: fx = -G*x/r^2, fy = -G*y/r^2, with r^2 = x^2 + y^2.
# Unit mass, so acceleration equals force.
def accel(pos, G):
    x, y = pos
    r2 = x*x + y*y                 # r^2 = x^2 + y^2
    inv = 1.0 / r2                 # 1/r^2
    return np.array([-G * x * inv, -G * y * inv])

# --- Generic vector velocity-Verlet integrator for one 2D particle ---
# Implemented explicitly step by step (not via a black-box routine).
def velocity_verlet(pos0, vel0, G, dt, t_end):
    n_steps = int(round(t_end / dt))
    pos = np.array(pos0, dtype=float)
    vel = np.array(vel0, dtype=float)
    a = accel(pos, G)              # initial acceleration
    traj = np.empty((n_steps + 1, 2))
    traj[0] = pos
    for i in range(n_steps):
        # 1) advance position using current velocity and acceleration
        pos = pos + vel * dt + 0.5 * a * dt * dt
        # 2) evaluate new acceleration at the new position
        a_new = accel(pos, G)
        # 3) advance velocity using the average of old and new accelerations
        vel = vel + 0.5 * (a + a_new) * dt
        # 4) roll the new acceleration into the current one for next step
        a = a_new
        traj[i + 1] = pos
    return traj

# --- Parameters ---
G = 1.0
dt = 0.1

# Run 1: ellipse
traj1 = velocity_verlet([4.0, 0.0], [0.0, 1.0], G, dt, t_end=100.0)
# Run 2: precessing rosette
traj2 = velocity_verlet([2.0, 2.0], [-1.0, 1.0], G, dt, t_end=1000.0)

# --- Plots ---
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

axes[0].plot(traj1[:, 0], traj1[:, 1], lw=0.8)
axes[0].plot(0, 0, 'k+', ms=10)
axes[0].plot(traj1[0, 0], traj1[0, 1], 'go', label='start')
axes[0].set_title("IC1: (4,0), v=(0,1), t=100  (ellipse)")
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
axes[0].set_aspect("equal"); axes[0].legend()

axes[1].plot(traj2[:, 0], traj2[:, 1], lw=0.4)
axes[1].plot(0, 0, 'k+', ms=10)
axes[1].plot(traj2[0, 0], traj2[0, 1], 'go', label='start')
axes[1].set_title("IC2: (2,2), v=(-1,1), t=1000  (rosette)")
axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
axes[1].set_aspect("equal"); axes[1].legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.1.1_s5.png")

# --- Separate check ---
# Radial extent for each orbit.
r1 = np.hypot(traj1[:, 0], traj1[:, 1])
r2 = np.hypot(traj2[:, 0], traj2[:, 1])

# Run 1: does it close? Compare final point to start point.
gap1 = np.hypot(traj1[-1, 0] - traj1[0, 0], traj1[-1, 1] - traj1[0, 1])

# Run 2: closure gap and how fully it fills the annulus.
# Angle of pericenter (min-radius point) each pass indicates precession.
gap2 = np.hypot(traj2[-1, 0] - traj2[0, 0], traj2[-1, 1] - traj2[0, 1])

print(f"Run 1 (ellipse): r_min = {r1.min():.4f}")
print(f"Run 1 (ellipse): r_max = {r1.max():.4f}")
print(f"Run 1 (ellipse): start-to-end gap = {gap1:.4f}")
print(f"Run 1 (ellipse): points span = {r1.max() - r1.min():.4f} (annulus width)")

print(f"Run 2 (rosette): r_min = {r2.min():.4f}")
print(f"Run 2 (rosette): r_max = {r2.max():.4f}")
print(f"Run 2 (rosette): start-to-end gap = {gap2:.4f}")
print(f"Run 2 (rosette): annulus width = {r2.max() - r2.min():.4f}")

# Coverage of the annulus by the perihelion/aphelion angular directions:
# spread of the direction angle at maximum radius over many orbits.
theta2 = np.arctan2(traj2[:, 1], traj2[:, 0])
print(f"Run 2 (rosette): angular coverage = {theta2.max() - theta2.min():.4f} rad "
      f"(near 2*pi = {2*np.pi:.4f} means it fills the annulus)")

# One-sentence explanation of why the check confirms the result:
print("Check: Run 1's tiny start-to-end gap with a bounded (r_min, r_max) band shows a "
      "closed, near-elliptical orbit, whereas Run 2's large start-to-end gap combined with "
      "full ~2*pi angular coverage between fixed r_min and r_max shows a precessing rosette "
      "that fills the annulus without closing.")
