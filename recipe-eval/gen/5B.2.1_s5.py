import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Physical parameters and initial conditions
# ----------------------------------------------------------------------
G  = 1.0          # coupling constant of the 1/r-type central force
dt = 0.01         # integration time step
T  = 100.0        # total integration time
nsteps = int(round(T / dt))

# State layout for both positions and velocities: [x1, y1, x2, y2]
q0 = np.array([2.0, 0.0, -2.0, 0.0])   # initial positions
v0 = np.array([0.0, 0.4,  0.0, -0.2])  # initial velocities (net y-momentum)


# ----------------------------------------------------------------------
# Acceleration = force (unit masses) for the two-body central force.
# f1 acts on body 1, f2 = -f1 acts on body 2 (equal and opposite).
# ----------------------------------------------------------------------
def accel(q):
    x1, y1, x2, y2 = q
    r12sq = (x1 - x2) ** 2 + (y1 - y2) ** 2   # squared separation
    # force on body 1 points from 1 toward 2
    a1x = G * (x2 - x1) / r12sq
    a1y = G * (y2 - y1) / r12sq
    # force on body 2 is the equal-and-opposite reaction
    a2x = G * (x1 - x2) / r12sq
    a2y = G * (y1 - y2) / r12sq
    return np.array([a1x, a1y, a2x, a2y])


# ----------------------------------------------------------------------
# Generic vector velocity-Verlet integrator, implemented explicitly.
# ----------------------------------------------------------------------
def velocity_verlet(q0, v0, dt, nsteps, accel):
    dim = len(q0)
    qs = np.empty((nsteps + 1, dim))
    vs = np.empty((nsteps + 1, dim))
    q = q0.copy()
    v = v0.copy()
    a = accel(q)                       # initial acceleration
    qs[0], vs[0] = q, v
    for i in range(nsteps):
        q = q + v * dt + 0.5 * a * dt ** 2   # position update
        a_new = accel(q)                     # new acceleration at updated position
        v = v + 0.5 * (a + a_new) * dt       # velocity update (average of old/new accel)
        a = a_new                            # carry acceleration forward
        qs[i + 1], vs[i + 1] = q, v
    return qs, vs


# ----------------------------------------------------------------------
# Run 1: raw initial velocities (nonzero total momentum -> drift)
# ----------------------------------------------------------------------
qs_raw, _ = velocity_verlet(q0, v0, dt, nsteps, accel)

# ----------------------------------------------------------------------
# Run 2: shift to the center-of-mass frame by subtracting the mean
# velocity per axis, so the total momentum is zero.
# ----------------------------------------------------------------------
vx_mean = (v0[0] + v0[2]) / 2.0   # mean x-velocity of the two bodies
vy_mean = (v0[1] + v0[3]) / 2.0   # mean y-velocity of the two bodies
v0_com = v0 - np.array([vx_mean, vy_mean, vx_mean, vy_mean])
qs_com, _ = velocity_verlet(q0, v0_com, dt, nsteps, accel)

# ----------------------------------------------------------------------
# Diagnostics
# ----------------------------------------------------------------------
px_raw = v0[0] + v0[2]
py_raw = v0[1] + v0[3]
px_com = v0_com[0] + v0_com[2]
py_com = v0_com[1] + v0_com[3]

print("Mean velocity subtracted per axis (vx_mean, vy_mean):", vx_mean, vy_mean)
print("Centered initial velocities [vx1, vy1, vx2, vy2]:", v0_com.tolist())
print("Total momentum RAW  (px, py):", px_raw, py_raw)
print("Total momentum COM  (px, py):", px_com, py_com)

# Center of mass over time (should drift in RAW, stay fixed in COM frame)
com_raw = np.column_stack(((qs_raw[:, 0] + qs_raw[:, 2]) / 2.0,
                           (qs_raw[:, 1] + qs_raw[:, 3]) / 2.0))
com_com = np.column_stack(((qs_com[:, 0] + qs_com[:, 2]) / 2.0,
                           (qs_com[:, 1] + qs_com[:, 3]) / 2.0))

print("RAW center of mass at t=0   (x, y):", com_raw[0, 0], com_raw[0, 1])
print("RAW center of mass at t=100 (x, y):", com_raw[-1, 0], com_raw[-1, 1])
print("COM-frame center of mass at t=0   (x, y):", com_com[0, 0], com_com[0, 1])
print("COM-frame center of mass at t=100 (x, y):", com_com[-1, 0], com_com[-1, 1])
print("RAW COM total displacement magnitude:",
      float(np.hypot(com_raw[-1, 0] - com_raw[0, 0], com_raw[-1, 1] - com_raw[0, 1])))
print("COM-frame COM total displacement magnitude:",
      float(np.hypot(com_com[-1, 0] - com_com[0, 0], com_com[-1, 1] - com_com[0, 1])))

# ----------------------------------------------------------------------
# Orbit plots: before (raw, drifting) and after (COM frame, stationary)
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(qs_raw[:, 0], qs_raw[:, 1], lw=0.8, label="body 1")
ax1.plot(qs_raw[:, 2], qs_raw[:, 3], lw=0.8, label="body 2")
ax1.plot(com_raw[:, 0], com_raw[:, 1], "k--", lw=0.8, label="center of mass")
ax1.set_title("Before COM shift (nonzero momentum: orbits drift)")
ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.axis("equal"); ax1.legend()

ax2.plot(qs_com[:, 0], qs_com[:, 1], lw=0.8, label="body 1")
ax2.plot(qs_com[:, 2], qs_com[:, 3], lw=0.8, label="body 2")
ax2.plot(com_com[:, 0], com_com[:, 1], "k.", ms=3, label="center of mass")
ax2.set_title("After COM shift (zero momentum: stationary pattern)")
ax2.set_xlabel("x"); ax2.set_ylabel("y"); ax2.axis("equal"); ax2.legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.2.1_s5.png")

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# The raw run's center of mass moves in a straight line while the bodies
# orbit it, whereas after subtracting the mean velocity the center of mass
# stays fixed and only the relative orbital motion remains, confirming that
# the drift was pure center-of-mass translation and the internal dynamics
# are unchanged by the frame shift.
print("Explanation: A nonzero total momentum makes the center of mass "
      "translate uniformly while the bodies orbit it; subtracting the mean "
      "velocity zeroes the momentum so the center of mass stays fixed and "
      "only the unchanged relative orbit remains, confirming the drift was "
      "solely center-of-mass motion.")
