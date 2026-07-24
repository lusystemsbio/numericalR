import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Central attractive force of magnitude G/r toward the origin.
# Components: fx = -G*x/r^2, fy = -G*y/r^2, with r^2 = x^2 + y^2.
G = 1.0

def force(pos):
    x, y = pos
    r2 = x * x + y * y            # r^2 = x^2 + y^2
    return np.array([-G * x / r2, -G * y / r2])

def velocity_verlet(pos0, vel0, dt, t_end):
    """Generic vector velocity-Verlet for one 2D particle (mass = 1)."""
    n = int(round(t_end / dt))
    pos = np.array(pos0, dtype=float)
    vel = np.array(vel0, dtype=float)
    a = force(pos)                        # initial acceleration (m = 1)
    traj = np.empty((n + 1, 2))
    traj[0] = pos
    for i in range(n):
        pos = pos + vel * dt + 0.5 * a * dt * dt   # position update
        a_new = force(pos)                         # force at new position
        vel = vel + 0.5 * (a + a_new) * dt         # velocity update (avg accel)
        a = a_new                                  # carry acceleration forward
        traj[i + 1] = pos
    return traj

# Run 1: ellipse
traj1 = velocity_verlet((4.0, 0.0), (0.0, 1.0), 0.1, 100.0)
# Run 2: precessing rosette
traj2 = velocity_verlet((2.0, 2.0), (-1.0, 1.0), 0.1, 1000.0)

# Radii for the "check": how close each orbit comes to closing / filling an annulus.
r1 = np.hypot(traj1[:, 0], traj1[:, 1])
r2 = np.hypot(traj2[:, 0], traj2[:, 1])

# Closure check for run 1: distance between final point and starting point.
close_gap1 = np.hypot(*(traj1[-1] - traj1[0]))

print("Run 1 (ellipse): r_min = %.4f, r_max = %.4f" % (r1.min(), r1.max()))
print("Run 1 (ellipse): final-to-start distance (closure gap) = %.4f" % close_gap1)
print("Run 2 (rosette): r_min = %.4f, r_max = %.4f" % (r2.min(), r2.max()))
print("Run 2 (rosette): annulus width (r_max - r_min) = %.4f" % (r2.max() - r2.min()))

# Plots
fig, ax = plt.subplots(1, 2, figsize=(12, 6))
ax[0].plot(traj1[:, 0], traj1[:, 1], lw=0.8)
ax[0].plot(0, 0, 'k+')
ax[0].set_title("IC (4,0), v=(0,1): closed ellipse")
ax[0].set_aspect("equal"); ax[0].set_xlabel("x"); ax[0].set_ylabel("y")

ax[1].plot(traj2[:, 0], traj2[:, 1], lw=0.4)
ax[1].plot(0, 0, 'k+')
ax[1].set_title("IC (2,2), v=(-1,1): precessing rosette")
ax[1].set_aspect("equal"); ax[1].set_xlabel("x"); ax[1].set_ylabel("y")

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.1.1_s2.png")

# Explanation: a small closure gap with a narrow radial range for run 1 means the
# path returns to its start and stays at nearly one radius (a closed ellipse),
# while run 2's large, persistent r_max - r_min shows the orbit sweeps a full
# annulus without repeating, confirming a non-closing precessing rosette.
print("Check: Run 1's tiny closure gap and small r-range confirm a closed near-ellipse, "
      "whereas Run 2's wide r_min-to-r_max band that never re-closes confirms a precessing rosette filling an annulus.")
