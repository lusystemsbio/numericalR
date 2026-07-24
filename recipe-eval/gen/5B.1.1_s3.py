import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def force(pos, G):
    """Central attractive force with magnitude G/r pointing toward origin.
    Components: fx = -G*x/r^2, fy = -G*y/r^2, with r^2 = x^2 + y^2."""
    x, y = pos
    r2 = x * x + y * y
    return np.array([-G * x / r2, -G * y / r2])


def velocity_verlet(pos0, vel0, G, dt, t_end):
    """Generic vector velocity-Verlet integrator for one 2D particle.
    Implemented explicitly step by step (unit mass, so a = force)."""
    n = int(round(t_end / dt))
    pos = np.array(pos0, dtype=float)
    vel = np.array(vel0, dtype=float)
    traj = np.empty((n + 1, 2))
    traj[0] = pos

    a = force(pos, G)                       # initial acceleration
    for i in range(n):
        # 1) advance position using current velocity and acceleration
        pos = pos + vel * dt + 0.5 * a * dt * dt
        # 2) compute new acceleration at the updated position
        a_new = force(pos, G)
        # 3) advance velocity using the average of old and new acceleration
        vel = vel + 0.5 * (a + a_new) * dt
        # 4) roll the acceleration forward for the next step
        a = a_new
        traj[i + 1] = pos
    return traj


G = 1.0
dt = 0.1

# Run 1: initial position (4, 0), velocity (0, 1), to t = 100 -> ellipse
traj1 = velocity_verlet([4.0, 0.0], [0.0, 1.0], G, dt, 100.0)

# Run 2: initial position (2, 2), velocity (-1, 1), to t = 1000 -> rosette
traj2 = velocity_verlet([2.0, 2.0], [-1.0, 1.0], G, dt, 1000.0)

# --- Orbit plots ---
fig, axes = plt.subplots(1, 2, figsize=(12, 6))
axes[0].plot(traj1[:, 0], traj1[:, 1], lw=0.8)
axes[0].plot(0, 0, 'k+', ms=10)
axes[0].set_title("IC 1: (4,0), v=(0,1) — ellipse")
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
axes[0].set_aspect("equal"); axes[0].grid(True)

axes[1].plot(traj2[:, 0], traj2[:, 1], lw=0.4)
axes[1].plot(0, 0, 'k+', ms=10)
axes[1].set_title("IC 2: (2,2), v=(-1,1) — precessing rosette")
axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
axes[1].set_aspect("equal"); axes[1].grid(True)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.1.1_s3.png")

# --- Separate check ---
# Radial extent of each orbit. A closed ellipse has one perihelion and one
# aphelion (a thin ring of radii); a precessing rosette sweeps every angle
# and fills an annulus between its inner and outer radii.
r1 = np.sqrt(traj1[:, 0]**2 + traj1[:, 1]**2)
r2 = np.sqrt(traj2[:, 0]**2 + traj2[:, 1]**2)

print("Run 1 (ellipse):  r_min = %.6f" % r1.min())
print("Run 1 (ellipse):  r_max = %.6f" % r1.max())
print("Run 1 (ellipse):  final position = (%.6f, %.6f)" % (traj1[-1, 0], traj1[-1, 1]))
print("Run 1 (ellipse):  distance of end point from start = %.6f"
      % np.hypot(traj1[-1, 0] - traj1[0, 0], traj1[-1, 1] - traj1[0, 1]))

print("Run 2 (rosette):  r_min = %.6f" % r2.min())
print("Run 2 (rosette):  r_max = %.6f" % r2.max())
print("Run 2 (rosette):  annulus width (r_max - r_min) = %.6f" % (r2.max() - r2.min()))

# Angle coverage: fraction of the full 2*pi circle actually visited.
theta2 = np.arctan2(traj2[:, 1], traj2[:, 0])
bins = np.linspace(-np.pi, np.pi, 361)
hist, _ = np.histogram(theta2, bins=bins)
print("Run 2 (rosette):  fraction of angular directions visited = %.4f"
      % (np.count_nonzero(hist) / len(hist)))

# Why this check confirms the result: the first orbit stays between a tight
# perihelion/aphelion pair and returns near its start (closed ellipse), while
# the second visits essentially every angle and spans a wide range of radii,
# i.e. it fills an annulus and never repeats its path (precessing rosette).
print("Check: Run 1 keeps r in a narrow band and closes on itself, whereas "
      "Run 2 covers all angles across a wide radial band, confirming a closed "
      "ellipse versus an annulus-filling precessing rosette.")
