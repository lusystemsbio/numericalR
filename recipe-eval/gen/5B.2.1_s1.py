import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Two equal (unit) masses attracting each other with a 1/r-type
# central force. State vector s = [x1, y1, x2, y2].
# Force on 1: f1 = G*(r2 - r1)/r12^2  (attractive, points 1 -> 2)
# Force on 2: equal and opposite.  Since m = 1, acceleration = force.
# ---------------------------------------------------------------

G = 1.0
dt = 0.01
t_end = 100.0
nsteps = int(round(t_end / dt))


def accel(s):
    """Return acceleration vector [a1x, a1y, a2x, a2y] for state s."""
    x1, y1, x2, y2 = s
    r12sq = (x1 - x2) ** 2 + (y1 - y2) ** 2  # squared separation
    # force on mass 1 (toward mass 2); mass 2 gets the opposite
    f1x = G * (x2 - x1) / r12sq
    f1y = G * (y2 - y1) / r12sq
    f2x = G * (x1 - x2) / r12sq
    f2y = G * (y1 - y2) / r12sq
    return np.array([f1x, f1y, f2x, f2y])


def verlet(s0, v0, nsteps, dt):
    """Generic vector velocity-Verlet integrator on all four coords."""
    S = np.zeros((nsteps + 1, 4))
    V = np.zeros((nsteps + 1, 4))
    S[0] = s0
    V[0] = v0
    a = accel(s0)                       # initial acceleration
    for i in range(nsteps):
        # half-step position update using current velocity and accel
        S[i + 1] = S[i] + V[i] * dt + 0.5 * a * dt * dt
        a_new = accel(S[i + 1])         # accel at the new position
        # velocity update using average of old and new accel
        V[i + 1] = V[i] + 0.5 * (a + a_new) * dt
        a = a_new                       # carry acceleration forward
    return S, V


def shift_to_com(S, V):
    """Subtract the mean velocity per axis (equal unit masses => CM frame)."""
    # mean velocity in x and in y from the two bodies (columns 0,2 and 1,3)
    Sc = S.copy()
    vx_mean = 0.5 * (V[:, 0] + V[:, 2])   # per-step mean x-velocity
    vy_mean = 0.5 * (V[:, 1] + V[:, 3])   # per-step mean y-velocity
    # remove the accumulated center-of-mass drift from each coordinate
    t = np.arange(S.shape[0]) * dt
    # constant drift velocity (mean over trajectory) times time
    dvx = vx_mean.mean()
    dvy = vy_mean.mean()
    Sc[:, 0] -= dvx * t
    Sc[:, 2] -= dvx * t
    Sc[:, 1] -= dvy * t
    Sc[:, 3] -= dvy * t
    return Sc, dvx, dvy


# --- initial conditions ---------------------------------------------------
s0 = np.array([2.0, 0.0, -2.0, 0.0])          # (x1, y1, x2, y2)
v0 = np.array([0.0, 0.4, 0.0, -0.2])          # (v1x, v1y, v2x, v2y)

# Run 1: keep the nonzero total momentum -> the pair drifts.
S1, V1 = verlet(s0, v0, nsteps, dt)

# total (mean) momentum per axis for run 1
px1 = 0.5 * (V1[0, 0] + V1[0, 2])
py1 = 0.5 * (V1[0, 1] + V1[0, 3])
print(f"Run 1 initial mean velocity (x, y): {px1:.6f}, {py1:.6f}")
print(f"Run 1 total momentum (x, y):        {V1[0,0]+V1[0,2]:.6f}, {V1[0,1]+V1[0,3]:.6f}")

# Shift run 1 into the center-of-mass frame (remove drift).
S1c, dvx1, dvy1 = shift_to_com(S1, V1)
print(f"Run 1 removed drift velocity (x, y): {dvx1:.6f}, {dvy1:.6f}")

# Run 2: center the velocities to zero total momentum before integrating.
vx_mean0 = 0.5 * (v0[0] + v0[2])
vy_mean0 = 0.5 * (v0[1] + v0[3])
v0_centered = v0.copy()
v0_centered[[0, 2]] -= vx_mean0    # remove mean x-velocity
v0_centered[[1, 3]] -= vy_mean0    # remove mean y-velocity
S2, V2 = verlet(s0, v0_centered, nsteps, dt)
print(f"Run 2 initial mean velocity (x, y): "
      f"{0.5*(v0_centered[0]+v0_centered[2]):.6f}, {0.5*(v0_centered[1]+v0_centered[3]):.6f}")
print(f"Run 2 total momentum (x, y):        "
      f"{v0_centered[0]+v0_centered[2]:.6f}, {v0_centered[1]+v0_centered[3]:.6f}")

# Center-of-mass displacement over the run (confirms drift vs. stationary).
comx1 = 0.5 * (S1[:, 0] + S1[:, 2])
comy1 = 0.5 * (S1[:, 1] + S1[:, 3])
comx2 = 0.5 * (S2[:, 0] + S2[:, 2])
comy2 = 0.5 * (S2[:, 1] + S2[:, 3])
print(f"Run 1 COM displacement to t=100 (x, y): "
      f"{comx1[-1]-comx1[0]:.6f}, {comy1[-1]-comy1[0]:.6f}")
print(f"Run 2 COM displacement to t=100 (x, y): "
      f"{comx2[-1]-comx2[0]:.6f}, {comy2[-1]-comy2[0]:.6f}")
comx1c = 0.5 * (S1c[:, 0] + S1c[:, 2])
comy1c = 0.5 * (S1c[:, 1] + S1c[:, 3])
print(f"Run 1 (shifted) COM displacement to t=100 (x, y): "
      f"{comx1c[-1]-comx1c[0]:.6f}, {comy1c[-1]-comy1c[0]:.6f}")

# --- plots ----------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15, 5))

ax[0].plot(S1[:, 0], S1[:, 1], label="body 1")
ax[0].plot(S1[:, 2], S1[:, 3], label="body 2")
ax[0].set_title("Run 1: before CM shift (drifting)")
ax[0].set_xlabel("x"); ax[0].set_ylabel("y")
ax[0].axis("equal"); ax[0].legend()

ax[1].plot(S1c[:, 0], S1c[:, 1], label="body 1")
ax[1].plot(S1c[:, 2], S1c[:, 3], label="body 2")
ax[1].set_title("Run 1: after CM shift (stationary)")
ax[1].set_xlabel("x"); ax[1].set_ylabel("y")
ax[1].axis("equal"); ax[1].legend()

ax[2].plot(S2[:, 0], S2[:, 1], label="body 1")
ax[2].plot(S2[:, 2], S2[:, 3], label="body 2")
ax[2].set_title("Run 2: zero-momentum init")
ax[2].set_xlabel("x"); ax[2].set_ylabel("y")
ax[2].axis("equal"); ax[2].legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.2.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the mutual central force is internal and cancels "
      "in pairs, the center of mass moves at constant velocity, so a nonzero "
      "total momentum makes the whole pattern translate steadily, and "
      "subtracting that mean velocity removes the drift to leave a stationary "
      "orbit about the common center of mass.")
