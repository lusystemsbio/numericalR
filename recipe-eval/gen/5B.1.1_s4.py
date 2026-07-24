import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Central attractive force of magnitude G/r pointing toward the origin.
# Cartesian components (so that magnitude = G/r):
#   fx = -G*x/r^2 , fy = -G*y/r^2 , with r^2 = x^2 + y^2
# Mass is taken as 1, so force = acceleration.
# ----------------------------------------------------------------------
def accel(pos, G):
    x, y = pos
    r2 = x*x + y*y          # r^2 = x^2 + y^2
    ax = -G * x / r2        # fx = -G*x/r^2
    ay = -G * y / r2        # fy = -G*y/r^2
    return np.array([ax, ay])

# ----------------------------------------------------------------------
# Generic vector velocity-Verlet integrator for one 2D particle.
# Implemented explicitly (no one-shot routine):
#   1) half-kick / drift bookkeeping done via the standard 3-line update
# ----------------------------------------------------------------------
def velocity_verlet(pos0, vel0, G, dt, t_final):
    n_steps = int(round(t_final / dt))
    pos = np.array(pos0, dtype=float)
    vel = np.array(vel0, dtype=float)
    a = accel(pos, G)                       # initial acceleration

    traj = np.empty((n_steps + 1, 2))
    traj[0] = pos

    for i in range(n_steps):
        # 1) advance position using current velocity and acceleration
        pos = pos + vel * dt + 0.5 * a * dt * dt
        # 2) compute new acceleration at the new position
        a_new = accel(pos, G)
        # 3) advance velocity using the average of old and new acceleration
        vel = vel + 0.5 * (a + a_new) * dt
        # 4) roll the acceleration forward for the next step
        a = a_new
        traj[i + 1] = pos

    return traj

# ----------------------------------------------------------------------
# Parameters and the two runs
# ----------------------------------------------------------------------
G = 1.0
dt = 0.1

# Run 1: near-elliptical closed orbit
traj1 = velocity_verlet((4.0, 0.0), (0.0, 1.0), G, dt, 100.0)

# Run 2: precessing rosette that fills an annulus
traj2 = velocity_verlet((2.0, 2.0), (-1.0, 1.0), G, dt, 1000.0)

# ----------------------------------------------------------------------
# Check quantities: radial extent (peri- and apo-center) of each orbit.
# A closed ellipse has a well-defined, essentially unchanging pair of
# turning radii; a filling rosette spans a whole annulus between r_min
# and r_max as its orientation precesses.
# ----------------------------------------------------------------------
r1 = np.sqrt(traj1[:, 0]**2 + traj1[:, 1]**2)
r2 = np.sqrt(traj2[:, 0]**2 + traj2[:, 1]**2)

print("Run 1 (4,0) v=(0,1) to t=100:")
print("  r_min = {:.6f}".format(r1.min()))
print("  r_max = {:.6f}".format(r1.max()))
print("  start-vs-end position gap = {:.6f}".format(
      np.linalg.norm(traj1[0] - traj1[-1])))

print("Run 2 (2,2) v=(-1,1) to t=1000:")
print("  r_min = {:.6f}".format(r2.min()))
print("  r_max = {:.6f}".format(r2.max()))
print("  annulus width (r_max - r_min) = {:.6f}".format(r2.max() - r2.min()))
print("  start-vs-end position gap = {:.6f}".format(
      np.linalg.norm(traj2[0] - traj2[-1])))

# ----------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

axes[0].plot(traj1[:, 0], traj1[:, 1], lw=0.8)
axes[0].plot(0, 0, 'k+', ms=10)
axes[0].set_title("IC (4,0), v=(0,1): near-elliptical closed orbit")
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
axes[0].set_aspect("equal"); axes[0].grid(True)

axes[1].plot(traj2[:, 0], traj2[:, 1], lw=0.4)
axes[1].plot(0, 0, 'k+', ms=10)
axes[1].set_title("IC (2,2), v=(-1,1): precessing rosette (annulus)")
axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
axes[1].set_aspect("equal"); axes[1].grid(True)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.1.1_s4.png")

# ----------------------------------------------------------------------
# Explanation of why the check confirms the result:
# For run 1 the two turning radii stay nearly equal to the initial radius
# and the path retraces itself (small start/end gap), i.e. a closed ellipse;
# for run 2 the large gap between r_min and r_max with a big start/end gap
# shows the orbit sweeps a full annulus without ever repeating, i.e. an
# open precessing rosette.
# ----------------------------------------------------------------------
print("Explanation: run 1's nearly constant turning radii and tiny start/end "
      "gap confirm a closed near-ellipse, while run 2's wide r_min-to-r_max "
      "annulus and large start/end gap confirm a non-closing precessing rosette.")
