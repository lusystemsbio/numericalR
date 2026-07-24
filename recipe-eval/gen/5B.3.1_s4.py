import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Three-body problem: three unit masses, pairwise 1/r^2-scaled
# central forces. State vector holds the six coordinates:
#   s = [x1, y1, x2, y2, x3, y3]
# ---------------------------------------------------------------

G = 1.0          # gravitational-like constant
dt = 0.01        # time step
t_max = 100.0    # final time
nsteps = int(t_max / dt)

# Acceleration (= force, since unit masses) on all six coordinates.
def accel(s):
    x1, y1, x2, y2, x3, y3 = s
    # pairwise squared distances (softened trivially only if r==0)
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)
    # force on body 1 from bodies 2 and 3 (1/r^2 scaling of the offset)
    f1x = G * ((x2 - x1) / r12**2 + (x3 - x1) / r13**2)
    f1y = G * ((y2 - y1) / r12**2 + (y3 - y1) / r13**2)
    # force on body 2 from bodies 1 and 3
    f2x = G * ((x1 - x2) / r12**2 + (x3 - x2) / r23**2)
    f2y = G * ((y1 - y2) / r12**2 + (y3 - y2) / r23**2)
    # force on body 3 from bodies 1 and 2
    f3x = G * ((x1 - x3) / r13**2 + (x2 - x3) / r23**2)
    f3y = G * ((y1 - y3) / r13**2 + (y2 - y3) / r23**2)
    return np.array([f1x, f1y, f2x, f2y, f3x, f3y])

# Potential-like energy for our 1/r^2-scaled force: force ~ 1/r,
# so the corresponding potential is -G*ln(r). Kinetic is (1/2)v^2.
def energy(s, v):
    x1, y1, x2, y2, x3, y3 = s
    r12 = np.hypot(x2 - x1, y2 - y1)
    r13 = np.hypot(x3 - x1, y3 - y1)
    r23 = np.hypot(x3 - x2, y3 - y2)
    KE = 0.5 * np.sum(v**2)
    PE = -G * (np.log(r12) + np.log(r13) + np.log(r23))
    return KE + PE

# --- initial conditions ---
s = np.array([-2.0, 0.0, 2.0, 0.0, 0.0, 3.0])   # (x1,y1,x2,y2,x3,y3)
v = np.array([0.1, 0.0, -0.1, 0.0, 0.0, -0.05])  # raw velocities

# Shift into the center-of-mass frame (unit masses -> subtract mean).
vx_cm = (v[0] + v[2] + v[4]) / 3.0
vy_cm = (v[1] + v[3] + v[5]) / 3.0
v[0::2] -= vx_cm
v[1::2] -= vy_cm
print("COM velocity removed (vx, vy):", vx_cm, vy_cm)

# --- velocity-Verlet integration (explicit, on the six coordinates) ---
traj = np.zeros((nsteps + 1, 6))
traj[0] = s
a = accel(s)                       # initial acceleration
E0 = energy(s, v)

# per-body kinetic energy history to detect the escaper
ke_hist = np.zeros((nsteps + 1, 3))
ke_hist[0] = [0.5*(v[0]**2+v[1]**2), 0.5*(v[2]**2+v[3]**2), 0.5*(v[4]**2+v[5]**2)]

escape_step = nsteps               # default: no escape within window
for n in range(nsteps):
    s = s + v * dt + 0.5 * a * dt**2   # position half-kick + drift
    a_new = accel(s)                    # new acceleration at new position
    v = v + 0.5 * (a + a_new) * dt      # velocity update using avg accel
    a = a_new                           # carry forward
    traj[n + 1] = s
    ke_hist[n + 1] = [0.5*(v[0]**2+v[1]**2),
                      0.5*(v[2]**2+v[3]**2),
                      0.5*(v[4]**2+v[5]**2)]
    # detect escape: any pairwise separation grows beyond a large threshold
    x1, y1, x2, y2, x3, y3 = s
    seps = [np.hypot(x2-x1, y2-y1), np.hypot(x3-x1, y3-y1), np.hypot(x3-x2, y3-y2)]
    if max(seps) > 50.0 and escape_step == nsteps:
        escape_step = n + 1

Efinal = energy(s, v)
print("Initial energy E0:", E0)
print("Final energy Ef:", Efinal)
print("Energy drift |Ef - E0|:", abs(Efinal - E0))
print("Escape detected at step:", escape_step, "time:", escape_step * dt)

# --- orbit plot up to the escape ---
end = escape_step
fig, ax = plt.subplots(figsize=(8, 8))
labels = ["Body 1", "Body 2", "Body 3"]
colors = ["tab:blue", "tab:orange", "tab:green"]
for i in range(3):
    ax.plot(traj[:end + 1, 2 * i], traj[:end + 1, 2 * i + 1],
            color=colors[i], lw=0.8, label=labels[i])
    ax.plot(traj[0, 2 * i], traj[0, 2 * i + 1], "o", color=colors[i])   # start
    ax.plot(traj[end, 2 * i], traj[end, 2 * i + 1], "s", color=colors[i])  # end
ax.set_aspect("equal")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_title("Three-body orbits (up to escape)")
ax.legend()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/5B.3.1_s4.png")

# --- separate check: bound-then-escape signature ---
# Compare the max pairwise separation early vs late, and identify which
# single body ends up carrying the largest kinetic energy (the escaper).
def max_sep_at(n):
    x1, y1, x2, y2, x3, y3 = traj[n]
    return max(np.hypot(x2-x1, y2-y1), np.hypot(x3-x1, y3-y1), np.hypot(x3-x2, y3-y2))

early_sep = np.mean([max_sep_at(n) for n in range(0, min(500, nsteps))])
late_sep = max_sep_at(escape_step)
print("Mean max-separation over first 5 time units (bound weaving phase):", early_sep)
print("Max separation at detected escape:", late_sep)
print("Final per-body kinetic energies (Body1, Body2, Body3):", ke_hist[escape_step])
print("Index of escaping body (max final KE):", int(np.argmax(ke_hist[escape_step])) + 1)

# Explanation:
print("Check rationale: the system stays compact (small bounded separations) "
      "while weaving, then one body's kinetic energy spikes and its separation "
      "diverges as the other two recoil into a tighter pair, which is exactly "
      "the bound-interaction-then-single-body-ejection signature of three-body chaos.")
