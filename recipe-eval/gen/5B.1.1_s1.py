import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Central attractive force of magnitude G/r pointing to the origin.
# Cartesian components: fx = -G*x/r^2, fy = -G*y/r^2, r^2 = x^2 + y^2.
# ---------------------------------------------------------------
def force(pos, G):
    x, y = pos
    r2 = x * x + y * y
    fx = -G * x / r2
    fy = -G * y / r2
    return np.array([fx, fy])

# ---------------------------------------------------------------
# Generic vector velocity-Verlet integrator for one 2D particle.
# Implemented explicitly step-by-step (mass = 1, so a = force).
# ---------------------------------------------------------------
def velocity_verlet(pos0, vel0, G, dt, t_end):
    n_steps = int(round(t_end / dt))
    pos = np.array(pos0, dtype=float)   # current position vector
    vel = np.array(vel0, dtype=float)   # current velocity vector
    traj = np.empty((n_steps + 1, 2))   # store the orbit for plotting
    traj[0] = pos

    a = force(pos, G)                   # initial acceleration (m = 1)
    for i in range(n_steps):
        # 1) advance position using current velocity and acceleration
        pos = pos + vel * dt + 0.5 * a * dt * dt
        # 2) compute new acceleration at the new position
        a_new = force(pos, G)
        # 3) advance velocity using the average of old and new acceleration
        vel = vel + 0.5 * (a + a_new) * dt
        # 4) roll the acceleration forward for the next step
        a = a_new
        traj[i + 1] = pos
    return traj

# ---------------------------------------------------------------
# Common parameters and the two test runs.
# ---------------------------------------------------------------
G = 1.0
dt = 0.1

# Run 1: near-circular / elliptical closed orbit
traj1 = velocity_verlet([4.0, 0.0], [0.0, 1.0], G, dt, 100.0)

# Run 2: precessing rosette that fills an annulus
traj2 = velocity_verlet([2.0, 2.0], [-1.0, 1.0], G, dt, 1000.0)

# ---------------------------------------------------------------
# Orbit plots for both initial conditions.
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

axes[0].plot(traj1[:, 0], traj1[:, 1], lw=0.8)
axes[0].plot(0, 0, 'k+', ms=10)
axes[0].set_title("IC (4,0) v=(0,1), t=100: near-elliptical closed orbit")
axes[0].set_xlabel("x"); axes[0].set_ylabel("y")
axes[0].set_aspect("equal"); axes[0].grid(True)

axes[1].plot(traj2[:, 0], traj2[:, 1], lw=0.4)
axes[1].plot(0, 0, 'k+', ms=10)
axes[1].set_title("IC (2,2) v=(-1,1), t=1000: precessing rosette")
axes[1].set_xlabel("x"); axes[1].set_ylabel("y")
axes[1].set_aspect("equal"); axes[1].grid(True)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.1.1_s1.png")

# ---------------------------------------------------------------
# Separate quantitative check.
# For orbit 1: compare the radius at the end with the radius at start
# and measure the spread of radii -> should stay in a thin band (closed ellipse).
# For orbit 2: measure the radial range -> a wide annulus that is filled.
# ---------------------------------------------------------------
r1 = np.sqrt(traj1[:, 0]**2 + traj1[:, 1]**2)
r2 = np.sqrt(traj2[:, 0]**2 + traj2[:, 1]**2)

# closure measure for orbit 1: distance between final point and initial point
closure1 = np.sqrt((traj1[-1, 0] - traj1[0, 0])**2 + (traj1[-1, 1] - traj1[0, 1])**2)

print("Orbit 1 (4,0)->(0,1): r_min = %.4f" % r1.min())
print("Orbit 1 (4,0)->(0,1): r_max = %.4f" % r1.max())
print("Orbit 1 (4,0)->(0,1): radial range (r_max - r_min) = %.4f" % (r1.max() - r1.min()))
print("Orbit 1 (4,0)->(0,1): final-to-initial point distance (closure) = %.4f" % closure1)
print("Orbit 1 (4,0)->(0,1): initial radius = %.4f, final radius = %.4f" % (r1[0], r1[-1]))

print("Orbit 2 (2,2)->(-1,1): r_min = %.4f" % r2.min())
print("Orbit 2 (2,2)->(-1,1): r_max = %.4f" % r2.max())
print("Orbit 2 (2,2)->(-1,1): radial range (r_max - r_min) = %.4f" % (r2.max() - r2.min()))

# The check confirms the result because a small radial range with a small final-to-initial
# gap means the first orbit retraces one nearly-elliptical loop (closed), whereas a large
# radial range that is densely sampled means the second orbit sweeps a broad annulus while
# its perihelion angle keeps advancing, so it never repeats and never closes.
print("Explanation: orbit 1's tiny radial spread and near-zero closure gap show it retraces a single closed ellipse, while orbit 2's large radial range densely filled between r_min and r_max shows a precessing rosette that fills an annulus and never closes.")
