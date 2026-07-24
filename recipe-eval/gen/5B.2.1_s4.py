import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Force model: two unit masses, 1/r-type central force ----
# f1x = G*(x2-x1)/r12^2, f2x = G*(x1-x2)/r12^2 (equal and opposite), same for y.
def accelerations(state, G):
    x1, y1, x2, y2 = state
    r12sq = (x1 - x2)**2 + (y1 - y2)**2   # squared separation
    f1x = G * (x2 - x1) / r12sq           # force on mass 1 (x)
    f1y = G * (y2 - y1) / r12sq           # force on mass 1 (y)
    f2x = G * (x1 - x2) / r12sq           # equal and opposite on mass 2 (x)
    f2y = G * (y1 - y2) / r12sq           # equal and opposite on mass 2 (y)
    # unit masses -> acceleration equals force
    return np.array([f1x, f1y, f2x, f2y])

# ---- Generic vector velocity-Verlet integrator on the 4 coordinates ----
def velocity_verlet(pos0, vel0, G, dt, nsteps):
    pos = np.array(pos0, dtype=float)
    vel = np.array(vel0, dtype=float)
    traj = np.empty((nsteps + 1, 4))
    traj[0] = pos
    acc = accelerations(pos, G)                # initial acceleration
    for i in range(nsteps):
        pos = pos + vel * dt + 0.5 * acc * dt**2   # position update
        acc_new = accelerations(pos, G)            # new acceleration at new position
        vel = vel + 0.5 * (acc + acc_new) * dt     # velocity update (average of accels)
        acc = acc_new                              # carry acceleration forward
        traj[i + 1] = pos
    return traj

# ---- Shift velocities to the center-of-mass frame (equal masses) ----
def shift_to_com_velocity(vel0):
    v = np.array(vel0, dtype=float)
    # mean velocity per axis over the two equal masses
    mean_vx = (v[0] + v[2]) / 2.0
    mean_vy = (v[1] + v[3]) / 2.0
    return np.array([v[0] - mean_vx, v[1] - mean_vy,
                     v[2] - mean_vx, v[3] - mean_vy])

# ---- Parameters ----
G = 1.0
dt = 0.01
t_final = 100.0
nsteps = int(round(t_final / dt))

pos0 = [2.0, 0.0, -2.0, 0.0]      # x1, y1, x2, y2
vel0 = [0.0, 0.4, 0.0, -0.2]      # vx1, vy1, vx2, vy2

# Run 1: original velocities (nonzero total momentum -> drift)
traj_drift = velocity_verlet(pos0, vel0, G, dt, nsteps)

# Run 2: velocities centered to zero total momentum (COM frame)
vel0_com = shift_to_com_velocity(vel0)
traj_com = velocity_verlet(pos0, vel0_com, G, dt, nsteps)

# ---- Diagnostics: total momentum (unit masses) ----
Px_orig = vel0[0] + vel0[2]
Py_orig = vel0[1] + vel0[3]
Px_com = vel0_com[0] + vel0_com[2]
Py_com = vel0_com[1] + vel0_com[3]
print(f"Original total momentum Px = {Px_orig}")
print(f"Original total momentum Py = {Py_orig}")
print(f"COM-frame velocities (vx1,vy1,vx2,vy2) = {tuple(vel0_com)}")
print(f"COM-frame total momentum Px = {Px_com}")
print(f"COM-frame total momentum Py = {Py_com}")

# ---- Center-of-mass drift check ----
# Common center of mass over time (equal masses): mean of the two positions per axis.
comx_drift = (traj_drift[:, 0] + traj_drift[:, 2]) / 2.0
comy_drift = (traj_drift[:, 1] + traj_drift[:, 3]) / 2.0
comx_com = (traj_com[:, 0] + traj_com[:, 2]) / 2.0
comy_com = (traj_com[:, 1] + traj_com[:, 3]) / 2.0

print(f"Drift run: COM start = ({comx_drift[0]:.6f}, {comy_drift[0]:.6f})")
print(f"Drift run: COM end   = ({comx_drift[-1]:.6f}, {comy_drift[-1]:.6f})")
print(f"Drift run: COM total displacement = "
      f"{np.hypot(comx_drift[-1]-comx_drift[0], comy_drift[-1]-comy_drift[0]):.6f}")
print(f"COM-frame run: COM start = ({comx_com[0]:.6f}, {comy_com[0]:.6f})")
print(f"COM-frame run: COM end   = ({comx_com[-1]:.6f}, {comy_com[-1]:.6f})")
print(f"COM-frame run: COM total displacement = "
      f"{np.hypot(comx_com[-1]-comx_com[0], comy_com[-1]-comy_com[0]):.6e}")

# ---- Orbit plots before and after shifting to the COM frame ----
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

ax = axes[0]
ax.plot(traj_drift[:, 0], traj_drift[:, 1], 'b-', lw=0.7, label="mass 1")
ax.plot(traj_drift[:, 2], traj_drift[:, 3], 'r-', lw=0.7, label="mass 2")
ax.plot(comx_drift, comy_drift, 'g--', lw=1.0, label="center of mass")
ax.set_title("Before COM shift (drifting)")
ax.set_xlabel("x"); ax.set_ylabel("y")
ax.axis("equal"); ax.legend()

ax = axes[1]
ax.plot(traj_com[:, 0], traj_com[:, 1], 'b-', lw=0.7, label="mass 1")
ax.plot(traj_com[:, 2], traj_com[:, 3], 'r-', lw=0.7, label="mass 2")
ax.plot(comx_com, comy_com, 'g--', lw=1.0, label="center of mass")
ax.set_title("After COM shift (stationary pattern)")
ax.set_xlabel("x"); ax.set_ylabel("y")
ax.axis("equal"); ax.legend()

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.2.1_s4.png")

# One-sentence explanation:
# The check confirms the result because with nonzero total momentum the center of
# mass moves in a straight line while the pair orbits (a drifting pattern), whereas
# subtracting the mean velocity zeroes the total momentum so the center of mass stays
# fixed and the orbit closes about it, proving the shift removed only the bulk drift.
print("Explanation: nonzero total momentum makes the COM translate steadily while the "
      "two bodies orbit, so removing the mean velocity zeroes the momentum and pins the "
      "COM in place, leaving a stationary orbital pattern and confirming only bulk drift was removed.")
