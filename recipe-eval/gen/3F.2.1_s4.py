import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y mutually repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Vector field f(X,Y) = (dX/dt, dY/dt)
def f(s):
    X, Y = s
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# Analytic Jacobian of f (used to find the stable direction at the saddle)
def jac(s):
    X, Y = s
    # d(dX)/dX, d(dX)/dY
    dxx = -kX
    hY = (Y / Yth) ** nY
    dxy = gX1 * (-(nY / Yth) * (Y / Yth) ** (nY - 1)) / (1.0 + hY) ** 2
    # d(dY)/dX, d(dY)/dY
    hX = (X / Xth) ** nX
    dyx = gY1 * (-(nX / Xth) * (X / Xth) ** (nX - 1)) / (1.0 + hX) ** 2
    dyy = -kY
    return np.array([[dxx, dxy], [dyx, dyy]])

# ----------------------------------------------------------------------
# 1. Locate the fixed points (two stable nodes + one saddle)
# ----------------------------------------------------------------------
saddle = fsolve(f, [190.0, 127.0], fprime=jac, full_output=False)
node_hiX = fsolve(f, [550.0, 40.0])   # X-dominant stable state
node_hiY = fsolve(f, [50.0, 350.0])   # Y-dominant stable state
print("Saddle point            :", saddle)
print("Stable node (X-dominant) :", node_hiX)
print("Stable node (Y-dominant) :", node_hiY)

# ----------------------------------------------------------------------
# 2. Find the STABLE eigen-direction at the saddle.
#    The separatrix is the saddle's stable manifold: forward flow moves
#    toward the saddle along it, so time-reversal moves outward along it.
# ----------------------------------------------------------------------
evals, evecs = np.linalg.eig(jac(saddle))
print("Jacobian eigenvalues at saddle:", evals)
stable_idx = np.argmin(evals.real)          # eigenvalue with most negative real part
v_stable = evecs[:, stable_idx].real
v_stable = v_stable / np.linalg.norm(v_stable)
print("Stable eigenvector (separatrix tangent):", v_stable)

# ----------------------------------------------------------------------
# 3. Integrate the TIME-REVERSED system  dX/dt = -f(X)  from ten points
#    seeded just beside the saddle along the stable eigenvector.
#    Under reversal the stable manifold becomes repelling, so the flow
#    carries these seeds outward and traces the separatrix.
# ----------------------------------------------------------------------
# Plotting / integration window
xmin, xmax, ymin, ymax = 0.0, 560.0, 0.0, 380.0

def rev_rhs(t, s):          # reversed field
    return -f(s)

# Stop integration when a trajectory leaves the window (avoids runaway)
def leave_box(t, s):
    X, Y = s
    inside = (X - xmin) * (xmax - X) * (Y - ymin) * (ymax - Y)
    return inside          # crosses zero when leaving the box
leave_box.terminal = True
leave_box.direction = -1

# Seed 5 offsets on each side of the saddle => 10 points near the saddle
distances = [0.5, 1.0, 1.5, 2.0, 2.5]
seeds = []
for d in distances:
    seeds.append(saddle + d * v_stable)   # + side of stable manifold
    seeds.append(saddle - d * v_stable)   # - side of stable manifold
print("Number of seed points near saddle:", len(seeds))

sep_branch_pos = []   # reversed trajectories going out one way
sep_branch_neg = []   # reversed trajectories going out the other way
for seed in seeds:
    sol = solve_ivp(rev_rhs, [0, 5000], seed, events=leave_box,
                    max_step=1.0, rtol=1e-9, atol=1e-9, dense_output=False)
    traj = sol.y.T
    # Classify by which side of the saddle it departed toward
    if np.dot(seed - saddle, v_stable) > 0:
        sep_branch_pos.append(traj)
    else:
        sep_branch_neg.append(traj)

# Build ONE continuous separatrix polyline: longest branch each side,
# joined through the saddle (neg branch reversed + saddle + pos branch).
def longest(branches):
    return max(branches, key=lambda a: len(a))
neg = longest(sep_branch_neg)
pos = longest(sep_branch_pos)
separatrix = np.vstack([neg[::-1], saddle[None, :], pos])
print("Separatrix polyline vertices:", separatrix.shape[0])

# ----------------------------------------------------------------------
# 4. Forward trajectories from random starts (sample dynamics)
# ----------------------------------------------------------------------
def fwd_rhs(t, s):
    return f(s)

rng = np.random.default_rng(0)
n_random = 40
fwd_trajs = []
attractor_label = []
for _ in range(n_random):
    s0 = np.array([rng.uniform(xmin, xmax), rng.uniform(ymin, ymax)])
    sol = solve_ivp(fwd_rhs, [0, 4000], s0, max_step=2.0, rtol=1e-8, atol=1e-8)
    traj = sol.y.T
    fwd_trajs.append(traj)
    endp = traj[-1]
    # which stable node did it reach?
    label = 0 if np.linalg.norm(endp - node_hiX) < np.linalg.norm(endp - node_hiY) else 1
    attractor_label.append(label)
print("Random forward starts reaching X-dominant node:", attractor_label.count(0))
print("Random forward starts reaching Y-dominant node:", attractor_label.count(1))

# ----------------------------------------------------------------------
# 5. CHECK: does any forward trajectory cross the separatrix?
#    Count segment-segment intersections between each forward trajectory
#    polyline and the separatrix polyline. A true basin boundary must
#    have ZERO crossings.
# ----------------------------------------------------------------------
def seg_intersect(p1, p2, p3, p4):
    # Do segment p1p2 and p3p4 properly intersect?
    def cross(o, a, b):
        return (a[0]-o[0])*(b[1]-o[1]) - (a[1]-o[1])*(b[0]-o[0])
    d1 = cross(p3, p4, p1)
    d2 = cross(p3, p4, p2)
    d3 = cross(p1, p2, p3)
    d4 = cross(p1, p2, p4)
    if ((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)):
        return True
    return False

total_crossings = 0
sep_seg = [(separatrix[i], separatrix[i+1]) for i in range(len(separatrix)-1)]
for traj in fwd_trajs:
    for i in range(len(traj)-1):
        a, b = traj[i], traj[i+1]
        for c, d in sep_seg:
            if seg_intersect(a, b, c, d):
                total_crossings += 1
print("Total forward-trajectory / separatrix crossings:", total_crossings)

# ----------------------------------------------------------------------
# 6. Phase-plane plot
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

# light vector-field arrows for context
xg, yg = np.meshgrid(np.linspace(xmin, xmax, 22), np.linspace(ymin, ymax, 22))
U = gX0 + gX1 / (1.0 + (yg / Yth) ** nY) - kX * xg
V = gY0 + gY1 / (1.0 + (xg / Xth) ** nX) - kY * yg
spd = np.hypot(U, V)
ax.quiver(xg, yg, U/spd, V/spd, color="0.8", pivot="mid", scale=40, width=0.002)

for k, traj in enumerate(fwd_trajs):
    col = "tab:blue" if attractor_label[k] == 0 else "tab:green"
    ax.plot(traj[:, 0], traj[:, 1], color=col, lw=0.8, alpha=0.6)

ax.plot(separatrix[:, 0], separatrix[:, 1], "r-", lw=2.5, label="separatrix (reversed flow)")
ax.plot(*node_hiX, "ko", ms=9, label="stable node")
ax.plot(*node_hiY, "ko", ms=9)
ax.plot(*saddle, "rs", ms=10, label="saddle")
ax.plot([s[0] for s in seeds], [s[1] for s in seeds], "r.", ms=4)

ax.set_xlim(xmin, xmax); ax.set_ylim(ymin, ymax)
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Toggle switch: separatrix via time-reversed flow")
ax.legend(loc="upper right", fontsize=8)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3F.2.1_s4.png")

# Why the check confirms the result:
print("Explanation: Because solution trajectories are unique and cannot cross, "
      "a curve that is itself a trajectory (the saddle's stable manifold) can "
      "never be crossed by any other trajectory; zero crossings therefore "
      "confirms it is the invariant boundary separating the two basins.")
