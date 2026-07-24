import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Two equal unit masses attracting via a 1/r-type central force:
#   f1x = G*(x2-x1)/r12^2 ,  f2x = G*(x1-x2)/r12^2   (equal & opposite)
#   likewise in y ; r12^2 = (x1-x2)^2 + (y1-y2)^2
# State vector q = [x1, y1, x2, y2]; masses are all 1 so a = F.
# We integrate with a generic vector velocity-Verlet scheme.
# ---------------------------------------------------------------

G = 1.0
dt = 0.01
t_end = 100.0
nsteps = int(round(t_end / dt))
m = np.array([1.0, 1.0, 1.0, 1.0])  # mass per coordinate (both bodies unit mass)

def forces(q):
    # Compute acceleration on each of the 4 coordinates.
    x1, y1, x2, y2 = q
    r12sq = (x1 - x2)**2 + (y1 - y2)**2   # squared separation
    fx = G * (x2 - x1) / r12sq            # force on body 1 in x
    fy = G * (y2 - y1) / r12sq            # force on body 1 in y
    # Body 2 feels the equal-and-opposite force.
    a = np.array([fx, fy, -fx, -fy]) / m
    return a

def velocity_verlet(q0, v0, dt, nsteps):
    # Generic vector velocity-Verlet on the full coordinate vector.
    q = np.array(q0, dtype=float)
    v = np.array(v0, dtype=float)
    traj = np.empty((nsteps + 1, 4))
    traj[0] = q
    a = forces(q)                         # initial acceleration
    for i in range(nsteps):
        v_half = v + 0.5 * dt * a         # half-step velocity kick
        q = q + dt * v_half               # full-step position drift
        a_new = forces(q)                 # recompute acceleration at new position
        v = v_half + 0.5 * dt * a_new     # complete the velocity step
        a = a_new
        traj[i + 1] = q
    return traj

def shift_to_com_frame(v0):
    # Subtract the mean velocity per axis so total momentum is zero.
    v = np.array(v0, dtype=float)
    vx_mean = (v[0] + v[2]) / 2.0         # mean x-velocity of the two bodies
    vy_mean = (v[1] + v[3]) / 2.0         # mean y-velocity of the two bodies
    v_shifted = v.copy()
    v_shifted[0] -= vx_mean
    v_shifted[2] -= vx_mean
    v_shifted[1] -= vy_mean
    v_shifted[3] -= vy_mean
    return v_shifted, vx_mean, vy_mean

# ---- Initial conditions -------------------------------------------------
q0 = [2.0, 0.0, -2.0, 0.0]
v0 = [0.0, 0.4, 0.0, -0.2]               # nonzero total momentum -> drift

# Run 1: as given (drifting).
traj_drift = velocity_verlet(q0, v0, dt, nsteps)

# Run 2: velocities centered to zero total momentum (COM frame).
v0_com, vx_mean, vy_mean = shift_to_com_frame(v0)
traj_com = velocity_verlet(q0, v0_com, dt, nsteps)

# ---- Numerical diagnostics ----------------------------------------------
px0 = v0[0] + v0[2]
py0 = v0[1] + v0[3]
print("Total px (original run):", px0)
print("Total py (original run):", py0)
print("Mean vx subtracted:", vx_mean)
print("Mean vy subtracted:", vy_mean)
px_com = v0_com[0] + v0_com[2]
py_com = v0_com[1] + v0_com[3]
print("Total px (COM run):", px_com)
print("Total py (COM run):", py_com)

# Center-of-mass position over time for each run.
com_drift = np.column_stack([(traj_drift[:, 0] + traj_drift[:, 2]) / 2.0,
                             (traj_drift[:, 1] + traj_drift[:, 3]) / 2.0])
com_com = np.column_stack([(traj_com[:, 0] + traj_com[:, 2]) / 2.0,
                           (traj_com[:, 1] + traj_com[:, 3]) / 2.0])

print("COM start (drift run): x=%.6f y=%.6f" % (com_drift[0, 0], com_drift[0, 1]))
print("COM end   (drift run): x=%.6f y=%.6f" % (com_drift[-1, 0], com_drift[-1, 1]))
print("COM displacement (drift run): dx=%.6f dy=%.6f"
      % (com_drift[-1, 0] - com_drift[0, 0], com_drift[-1, 1] - com_drift[0, 1]))
print("COM start (COM run): x=%.6f y=%.6f" % (com_com[0, 0], com_com[0, 1]))
print("COM end   (COM run): x=%.6f y=%.6f" % (com_com[-1, 0], com_com[-1, 1]))
print("COM max drift from origin (COM run): %.6e"
      % np.max(np.hypot(com_com[:, 0], com_com[:, 1])))

# ---- Plots: before vs after shifting to the COM frame -------------------
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

ax = axes[0]
ax.plot(traj_drift[:, 0], traj_drift[:, 1], lw=0.8, label="body 1")
ax.plot(traj_drift[:, 2], traj_drift[:, 3], lw=0.8, label="body 2")
ax.plot(com_drift[:, 0], com_drift[:, 1], "k--", lw=1.0, label="center of mass")
ax.set_title("Before: nonzero momentum (orbit + drift)")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.axis("equal"); ax.legend()

ax = axes[1]
ax.plot(traj_com[:, 0], traj_com[:, 1], lw=0.8, label="body 1")
ax.plot(traj_com[:, 2], traj_com[:, 3], lw=0.8, label="body 2")
ax.plot(com_com[:, 0], com_com[:, 1], "k+", ms=8, label="center of mass")
ax.set_title("After: COM frame (stationary pattern)")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.axis("equal"); ax.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.2.1_s3.png")

# One-sentence explanation of why this check confirms the result:
print("Explanation: The drifting run shows the orbit pattern translating across the "
      "plane (COM moves), while subtracting the mean velocity zeroes the total momentum "
      "so the COM stays fixed and the same orbit becomes a stationary closed pattern, "
      "confirming that the shift only removes bulk translation and leaves the relative motion intact.")
