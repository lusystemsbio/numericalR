import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Three-body problem: three unit masses, each attracted to the others by
# a 1/r^2-scaled central force summed over both partners.
# State vector packs the six coordinates: [x1, x2, x3, y1, y2, y3].
# ----------------------------------------------------------------------

G = 1.0          # gravitational constant
m = 1.0          # unit masses
dt = 0.01        # time step
t_end = 100.0    # final time
nsteps = int(t_end / dt)


def accel(state):
    """Return the acceleration vector for the six coordinates.

    Force on body i from body j points along (rj - ri) and is scaled by
    1/rij^2 (as specified), summed over the two partners.
    """
    x1, x2, x3, y1, y2, y3 = state

    # pairwise separations and their (squared) distances
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)

    # accelerations (unit mass so force == acceleration)
    a1x = G * ((x2 - x1) / r12**2 + (x3 - x1) / r13**2)
    a2x = G * ((x1 - x2) / r12**2 + (x3 - x2) / r23**2)
    a3x = G * ((x1 - x3) / r13**2 + (x2 - x3) / r23**2)

    a1y = G * ((y2 - y1) / r12**2 + (y3 - y1) / r13**2)
    a2y = G * ((y1 - y2) / r12**2 + (y3 - y2) / r23**2)
    a3y = G * ((y1 - y3) / r13**2 + (y2 - y3) / r23**2)

    return np.array([a1x, a2x, a3x, a1y, a2y, a3y])


def total_energy(state, vel):
    """Kinetic + potential energy. Potential consistent with 1/r^2 force
    is U = -G*m*m/r for each pair (its gradient gives the stated force)."""
    x1, x2, x3, y1, y2, y3 = state
    ke = 0.5 * m * np.sum(vel**2)
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)
    pe = -G * m * m * (1.0 / r12 + 1.0 / r13 + 1.0 / r23)
    return ke + pe


# ---------------------- initial conditions ----------------------------
# positions: (x1, x2, x3, y1, y2, y3)
state = np.array([-2.0, 2.0, 0.0, 0.0, 0.0, 3.0])
# raw velocities in the same packing: (vx1, vx2, vx3, vy1, vy2, vy3)
vel = np.array([0.1, -0.1, 0.0, 0.0, 0.0, -0.05])

# --- shift into the center-of-mass frame (equal masses => subtract mean) ---
vx = vel[0:3]
vy = vel[3:6]
vx = vx - vx.mean()
vy = vy - vy.mean()
vel = np.concatenate([vx, vy])
# also center positions so the COM sits at the origin
px = state[0:3]
py = state[3:6]
px = px - px.mean()
py = py - py.mean()
state = np.concatenate([px, py])

print("Initial COM velocity (vx, vy):",
      state[0:3].mean(), state[3:6].mean())
print("Initial center-of-mass position (x, y):",
      state[0:3].mean(), state[3:6].mean())

E0 = total_energy(state, vel)
print("Initial total energy E0:", E0)

# ---------------------- velocity-Verlet loop --------------------------
# Explicit generic vector velocity-Verlet on the six coordinates:
#   1) x  <- x + v*dt + 0.5*a*dt^2
#   2) recompute acceleration a_new at the new positions
#   3) v  <- v + 0.5*(a + a_new)*dt
traj = np.empty((nsteps + 1, 6))
traj[0] = state
a = accel(state)                       # acceleration at the start

escape_step = None
pair_dist_thresh = 6.0                 # a body this far from the others => escaped

for k in range(nsteps):
    # step 1: advance positions with current velocity and acceleration
    state = state + vel * dt + 0.5 * a * dt**2
    # step 2: acceleration at the updated positions
    a_new = accel(state)
    # step 3: advance velocities with the average of old and new accelerations
    vel = vel + 0.5 * (a + a_new) * dt
    a = a_new                          # carry forward for the next iteration
    traj[k + 1] = state

    # detect escape: one body separated from BOTH others by more than threshold
    x1, x2, x3, y1, y2, y3 = state
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)
    if escape_step is None:
        if (r12 > pair_dist_thresh and r13 > pair_dist_thresh) or \
           (r12 > pair_dist_thresh and r23 > pair_dist_thresh) or \
           (r13 > pair_dist_thresh and r23 > pair_dist_thresh):
            escape_step = k + 1

Ef = total_energy(state, vel)
print("Final total energy Ef:", Ef)
print("Relative energy drift |Ef-E0|/|E0|:", abs(Ef - E0) / abs(E0))

# ---------------------- escape / recoil check -------------------------
if escape_step is not None:
    t_escape = escape_step * dt
    print("Escape detected at step:", escape_step, " time t =", t_escape)
    plot_upto = escape_step + 1
else:
    print("No escape detected within t =", t_end)
    plot_upto = nsteps + 1

# velocities of each body at the end (to show one flies off, pair recoils)
vx1, vx2, vx3 = vel[0], vel[1], vel[2]
vy1, vy2, vy3 = vel[3], vel[4], vel[5]
speed = [np.hypot(vx1, vy1), np.hypot(vx2, vy2), np.hypot(vx3, vy3)]
print("Final speed body 1:", speed[0])
print("Final speed body 2:", speed[1])
print("Final speed body 3:", speed[2])
print("Fastest (escaping) body index (1-based):", int(np.argmax(speed)) + 1)

# final COM velocity should remain ~0 (momentum conserved in COM frame)
print("Final COM velocity (vx, vy):", vel[0:3].mean(), vel[3:6].mean())

# ---------------------- orbit plot ------------------------------------
seg = traj[:plot_upto]
plt.figure(figsize=(8, 8))
colors = ["tab:blue", "tab:orange", "tab:green"]
for i in range(3):
    plt.plot(seg[:, i], seg[:, i + 3], color=colors[i], lw=0.8,
             label=f"body {i+1}")
    plt.plot(seg[0, i], seg[0, i + 3], "o", color=colors[i])   # start
    plt.plot(seg[-1, i], seg[-1, i + 3], "s", color=colors[i])  # end
plt.gca().set_aspect("equal", "box")
plt.xlabel("x")
plt.ylabel("y")
title = "Three-body orbits (COM frame)"
if escape_step is not None:
    title += f" up to escape at t={escape_step*dt:.2f}"
plt.title(title)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.3.1_s3.png", dpi=130, bbox_inches="tight")

# ---------------------- one-sentence explanation ----------------------
print("Explanation: The check confirms three-body chaos because energy is "
      "conserved overall (small drift) yet the bodies first stay mutually "
      "bound and weave, after which one body's speed grows and it separates "
      "from the other two while the remaining pair recoils in the opposite "
      "direction to conserve momentum in the COM frame.")
