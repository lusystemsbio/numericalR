import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Force model: two unit masses, 1/r-type central attraction -----
# state vector q = [x1, y1, x2, y2]; masses are all 1 so a = F.
def accelerations(q, G):
    x1, y1, x2, y2 = q
    r12sq = (x1 - x2)**2 + (y1 - y2)**2   # squared separation
    # force on 1 points toward 2, force on 2 is equal and opposite
    f1x = G * (x2 - x1) / r12sq
    f1y = G * (y2 - y1) / r12sq
    f2x = G * (x1 - x2) / r12sq
    f2y = G * (y1 - y2) / r12sq
    return np.array([f1x, f1y, f2x, f2y])

# ----- Generic vector velocity-Verlet integrator on the 4 coordinates -----
def velocity_verlet(q0, v0, G, dt, tmax):
    nsteps = int(round(tmax / dt))
    q = np.array(q0, dtype=float)
    v = np.array(v0, dtype=float)
    traj = np.zeros((nsteps + 1, 4))
    traj[0] = q
    a = accelerations(q, G)          # initial acceleration
    for n in range(nsteps):
        q = q + v * dt + 0.5 * a * dt**2   # position half-kick-drift update
        a_new = accelerations(q, G)         # new acceleration at new position
        v = v + 0.5 * (a + a_new) * dt      # velocity update uses avg accel
        a = a_new                           # carry acceleration forward
        traj[n + 1] = q
    return traj

# ----- Shift to center-of-mass frame: subtract mean velocity per axis -----
# Equal masses => mean over the two bodies is the COM velocity component.
def center_of_mass_shift(v):
    v = np.array(v, dtype=float)
    vx_mean = (v[0] + v[2]) / 2.0   # mean x-velocity (COM x-velocity)
    vy_mean = (v[1] + v[3]) / 2.0   # mean y-velocity (COM y-velocity)
    return np.array([v[0] - vx_mean, v[1] - vy_mean,
                     v[2] - vx_mean, v[3] - vy_mean])

# ----- Parameters and initial conditions -----
G = 1.0
dt = 0.01
tmax = 100.0
q0 = [2.0, 0.0, -2.0, 0.0]
v0 = [0.0, 0.4, 0.0, -0.2]

# Run 1: original velocities (nonzero total momentum -> drift)
traj_before = velocity_verlet(q0, v0, G, dt, tmax)

# Run 2: velocities centered to zero total momentum (COM frame)
v0_shifted = center_of_mass_shift(v0)
traj_after = velocity_verlet(q0, v0_shifted, G, dt, tmax)

# ----- Report numerical results -----
p_before = np.array([v0[0] + v0[2], v0[1] + v0[3]])          # total momentum (masses=1)
p_after = np.array([v0_shifted[0] + v0_shifted[2],
                    v0_shifted[1] + v0_shifted[3]])
print("Original velocities [v1x, v1y, v2x, v2y]:", v0)
print("Total momentum before shift (px, py):", p_before[0], p_before[1])
print("Shifted velocities [v1x, v1y, v2x, v2y]:", list(v0_shifted))
print("Total momentum after shift (px, py):", p_after[0], p_after[1])

# COM position over time (mean of the two bodies), for each run
com_before = np.column_stack([(traj_before[:, 0] + traj_before[:, 2]) / 2.0,
                              (traj_before[:, 1] + traj_before[:, 3]) / 2.0])
com_after = np.column_stack([(traj_after[:, 0] + traj_after[:, 2]) / 2.0,
                             (traj_after[:, 1] + traj_after[:, 3]) / 2.0])
print("COM initial (before shift):", com_before[0, 0], com_before[0, 1])
print("COM final   (before shift):", com_before[-1, 0], com_before[-1, 1])
print("COM displacement (before shift):",
      com_before[-1, 0] - com_before[0, 0], com_before[-1, 1] - com_before[0, 1])
print("COM initial (after shift):", com_after[0, 0], com_after[0, 1])
print("COM final   (after shift):", com_after[-1, 0], com_after[-1, 1])
print("COM displacement (after shift):",
      com_after[-1, 0] - com_after[0, 0], com_after[-1, 1] - com_after[0, 1])
print("Max |COM drift| before shift:",
      np.max(np.abs(com_before - com_before[0])))
print("Max |COM drift| after shift:",
      np.max(np.abs(com_after - com_after[0])))

# ----- Orbit plots: before vs after COM-frame shift -----
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

ax1.plot(traj_before[:, 0], traj_before[:, 1], 'b-', lw=0.7, label="mass 1")
ax1.plot(traj_before[:, 2], traj_before[:, 3], 'r-', lw=0.7, label="mass 2")
ax1.plot(com_before[:, 0], com_before[:, 1], 'k--', lw=0.9, label="COM")
ax1.set_title("Before COM shift (nonzero momentum: drifts)")
ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.axis("equal"); ax1.legend()

ax2.plot(traj_after[:, 0], traj_after[:, 1], 'b-', lw=0.7, label="mass 1")
ax2.plot(traj_after[:, 2], traj_after[:, 3], 'r-', lw=0.7, label="mass 2")
ax2.plot(com_after[:, 0], com_after[:, 1], 'k.', ms=2, label="COM")
ax2.set_title("After COM shift (zero momentum: stationary pattern)")
ax2.set_xlabel("x"); ax2.set_ylabel("y"); ax2.axis("equal"); ax2.legend()

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.2.1_s2.png")

# Explanation:
# The check confirms the result because a nonzero total momentum makes the COM
# translate uniformly (large COM drift), so the orbit smears across the plane,
# whereas subtracting the mean velocity zeroes total momentum and pins the COM
# (near-zero drift), leaving a closed, stationary orbit pattern about it.
print("Explanation: Zero total momentum pins the center of mass in place, so the "
      "smeared drifting orbit collapses to a stationary closed pattern about the COM.")
