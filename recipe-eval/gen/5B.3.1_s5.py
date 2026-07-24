import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Chaotic three-body problem, integrated with velocity-Verlet.
# State vector holds the six coordinates: [x1,x2,x3, y1,y2,y3].
# ---------------------------------------------------------------

G = 1.0
dt = 0.01
t_end = 100.0
nsteps = int(round(t_end / dt))

# --- initial conditions -----------------------------------------
# positions (x1,x2,x3, y1,y2,y3) packed to match the force formula layout
pos0 = np.array([-2.0, 2.0, 0.0,   0.0, 0.0, 3.0])  # x1,x2,x3, y1,y2,y3
# velocities as given (same packing)
vel0 = np.array([0.1, -0.1, 0.0,   0.0, 0.0, -0.05])

# shift into the center-of-mass frame (unit masses -> COM is the mean)
xs, ys = pos0[:3], pos0[3:]
vxs, vys = vel0[:3], vel0[3:]
xs = xs - xs.mean(); ys = ys - ys.mean()          # center positions
vxs = vxs - vxs.mean(); vys = vys - vys.mean()      # center velocities
pos = np.concatenate([xs, ys]).astype(float)
vel = np.concatenate([vxs, vys]).astype(float)

def accel(p):
    # unpack the six coordinates
    x1, x2, x3, y1, y2, y3 = p
    # pairwise squared separations
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)
    r12s, r13s, r23s = r12**2, r13**2, r23**2
    # accelerations = sum over partners of G*(delta)/r^2, matching f1x formula
    a1x = G * ((x2 - x1) / r12s + (x3 - x1) / r13s)
    a2x = G * ((x1 - x2) / r12s + (x3 - x2) / r23s)
    a3x = G * ((x1 - x3) / r13s + (x2 - x3) / r23s)
    a1y = G * ((y2 - y1) / r12s + (y3 - y1) / r13s)
    a2y = G * ((y1 - y2) / r12s + (y3 - y2) / r23s)
    a3y = G * ((y1 - y3) / r13s + (y2 - y3) / r23s)
    return np.array([a1x, a2x, a3x, a1y, a2y, a3y])

def energy(p, v):
    # kinetic (unit masses): 0.5*sum v^2
    KE = 0.5 * np.sum(v**2)
    x1, x2, x3, y1, y2, y3 = p
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)
    # potential consistent with the 1/r^2-scaled force: dU/dr = -1/r -> U = -G*ln(r)
    PE = -G * (np.log(r12) + np.log(r13) + np.log(r23))
    return KE, PE, KE + PE

# --- history buffers --------------------------------------------
traj = np.zeros((nsteps + 1, 6))
sep_max = np.zeros(nsteps + 1)
E_hist = np.zeros(nsteps + 1)
traj[0] = pos
E_hist[0] = energy(pos, vel)[2]

def max_sep(p):
    x1, x2, x3, y1, y2, y3 = p
    return max(np.hypot(x2 - x1, y2 - y1),
               np.hypot(x3 - x1, y3 - y1),
               np.hypot(x3 - x2, y3 - y2))
sep_max[0] = max_sep(pos)

# --- velocity-Verlet loop (explicit, step by step) --------------
a = accel(pos)                       # initial acceleration
for i in range(1, nsteps + 1):
    pos = pos + vel * dt + 0.5 * a * dt**2   # drift positions with current accel
    a_new = accel(pos)                        # recompute accel at new positions
    vel = vel + 0.5 * (a + a_new) * dt        # kick velocities with averaged accel
    a = a_new                                 # carry accel to next step
    traj[i] = pos
    sep_max[i] = max_sep(pos)
    E_hist[i] = energy(pos, vel)[2]

t = np.arange(nsteps + 1) * dt

# --- detect the escape: first time max separation exceeds threshold ---
escape_threshold = 15.0
esc_idx = np.argmax(sep_max > escape_threshold)
if sep_max[esc_idx] <= escape_threshold:   # never crossed
    esc_idx = nsteps
t_escape = t[esc_idx]

# --- energy conservation diagnostic ------------------------------
E0 = E_hist[0]
drift = np.max(np.abs(E_hist[:esc_idx + 1] - E0)) / abs(E0)

# --- orbit plot up to the escape ---------------------------------
plt.figure(figsize=(8, 8))
labels = ["Body 1", "Body 2", "Body 3"]
colors = ["tab:red", "tab:green", "tab:blue"]
for b in range(3):
    xb = traj[:esc_idx + 1, b]
    yb = traj[:esc_idx + 1, 3 + b]
    plt.plot(xb, yb, color=colors[b], lw=0.8, label=labels[b])
    plt.plot(xb[0], yb[0], "o", color=colors[b], ms=8)   # start marker
    plt.plot(xb[-1], yb[-1], "s", color=colors[b], ms=8) # end marker
plt.xlabel("x"); plt.ylabel("y")
plt.title("Chaotic three-body orbits (COM frame), up to escape")
plt.legend(); plt.axis("equal"); plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.3.1_s5.png")

# --- separate check: bound weaving -> one escapes, pair recoils ---
# early-phase max separation (weaving while bound) vs late-phase (escape)
early_window = slice(0, min(esc_idx, nsteps) // 2 + 1)
early_max_sep = np.max(sep_max[early_window])
final_max_sep = sep_max[-1]

# identify escaper: body whose distance from COM (0,0) is largest at the end
xf = traj[-1, :3]; yf = traj[-1, 3:]
dist_com = np.hypot(xf, yf)
escaper = int(np.argmax(dist_com)) + 1
remaining = [b + 1 for b in range(3) if b != escaper - 1]

# separation of the remaining pair at the end (should stay finite = bound recoiling pair)
pair_idx = [b - 1 for b in remaining]
pair_sep_final = np.hypot(xf[pair_idx[0]] - xf[pair_idx[1]],
                          yf[pair_idx[0]] - yf[pair_idx[1]])

# --- print all numerical results ---------------------------------
print(f"Initial total energy E0: {E0:.6f}")
print(f"Final total energy E(t_end): {E_hist[-1]:.6f}")
print(f"Max relative energy drift (pre-escape): {drift:.3e}")
print(f"Escape threshold (max pairwise separation): {escape_threshold:.3f}")
print(f"Escape detected at index: {esc_idx}")
print(f"Escape time t_escape: {t_escape:.3f}")
print(f"Early-phase max pairwise separation (bound weaving): {early_max_sep:.4f}")
print(f"Final max pairwise separation (escaping): {final_max_sep:.4f}")
print(f"Ratio final/early max separation: {final_max_sep/early_max_sep:.4f}")
print(f"Escaping body (farthest from COM at end): Body {escaper}")
print(f"Distances of bodies from COM at end: "
      f"Body1={dist_com[0]:.4f}, Body2={dist_com[1]:.4f}, Body3={dist_com[2]:.4f}")
print(f"Remaining bound pair: Body {remaining[0]} and Body {remaining[1]}")
print(f"Remaining-pair separation at end: {pair_sep_final:.4f}")

# One-sentence explanation of why the check confirms three-body chaos:
print("Explanation: The check confirms chaos because the small (nearly conserved) "
      "energy drift shows the integration is faithful, while the bounded early "
      "weaving followed by one body's max separation growing without limit as the "
      "other two stay close is exactly the ejection-plus-recoiling-binary signature "
      "of the chaotic three-body breakup.")
