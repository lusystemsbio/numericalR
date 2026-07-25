import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# ----------------------------------------------------------------------
# 6D.1 model: genetic toggle switch SDE (two mutually repressing genes).
#   dX = ( b / (1 + Y^n) - X ) dt + sigma dW_x
#   dY = ( b / (1 + X^n) - Y ) dt + sigma dW_y
# Two stable states: (high X, low Y) and (low X, high Y).
# Same trajectory as 6D.1: seed 3, b = 20, dt = 0.01, integrate to t = 1000.
# ----------------------------------------------------------------------
rng   = np.random.default_rng(3)      # seed 3
b     = 20.0                          # max production rate
n     = 4.0                           # Hill coefficient
sigma = 2.0                           # noise amplitude
dt    = 0.01
T     = 1000.0
nsteps = int(round(T / dt))
sqdt  = np.sqrt(dt)

# Euler-Maruyama integration of the toggle-switch SDE
X = np.empty(nsteps + 1)
Y = np.empty(nsteps + 1)
X[0], Y[0] = b, 0.0                    # start in the "high X" state
for k in range(nsteps):
    fx = b / (1.0 + Y[k]**n) - X[k]
    fy = b / (1.0 + X[k]**n) - Y[k]
    X[k+1] = X[k] + fx*dt + sigma*sqdt*rng.standard_normal()
    Y[k+1] = Y[k] + fy*dt + sigma*sqdt*rng.standard_normal()
X = np.clip(X, 0, None)               # concentrations stay non-negative
Y = np.clip(Y, 0, None)

pts = np.column_stack([X, Y])

# ----------------------------------------------------------------------
# Split the plane along the symmetry line Y = X to label the two clouds:
#   state 0 : below the line (Y < X)  -> "high X"
#   state 1 : above the line (Y > X)  -> "high Y"
# These labels are used ONLY to fit each state's ellipse from its cloud.
# ----------------------------------------------------------------------
below = Y < X                         # rough membership for fitting
clouds = {0: pts[below], 1: pts[~below]}

ENLARGE = 1.5                         # enlarge each ellipse's axes by 1.5x

ellipses = {}                         # state -> (mean, inv_cov_scaled)
ell_patches = []
for s, cloud in clouds.items():
    mu  = cloud.mean(axis=0)                       # ellipse center
    cov = np.cov(cloud, rowvar=False)              # 2x2 covariance
    # eigen-decomposition of covariance gives ellipse axes & orientation
    evals, evecs = np.linalg.eigh(cov)             # evals = variances along axes
    # semi-axis lengths (1 std) enlarged by the factor 1.5
    axis_len = ENLARGE * np.sqrt(evals)
    # A point p is inside the ellipse iff, in the eigenbasis,
    #   sum_i ( (p-mu).v_i / axis_len_i )^2 <= 1
    # Build the scaled inverse-shape matrix M so test is (p-mu) M (p-mu) <= 1:
    #   M = V diag(1/axis_len^2) V^T
    M = evecs @ np.diag(1.0 / axis_len**2) @ evecs.T
    ellipses[s] = (mu, M)
    # store a drawable patch (angle & full-axis lengths for matplotlib)
    order = np.argsort(evals)[::-1]
    major, minor = axis_len[order]
    angle = np.degrees(np.arctan2(evecs[1, order[0]], evecs[0, order[0]]))
    ell_patches.append((mu, 2*major, 2*minor, angle, s))

def inside(p, s):
    """Return True if point p lies inside state s's enlarged ellipse."""
    mu, M = ellipses[s]
    d = p - mu
    return d @ M @ d <= 1.0

# ----------------------------------------------------------------------
# Assign each point to a state:
#   - if inside its own-side ellipse -> that state
#   - otherwise keep the previous assignment (point is between clouds)
# A transition is a change in the assigned state between consecutive steps.
# ----------------------------------------------------------------------
assigned = np.empty(nsteps + 1, dtype=int)
# initialize from side of Y=X line
cur = 0 if pts[0, 1] < pts[0, 0] else 1
for k in range(nsteps + 1):
    p = pts[k]
    # side of the Y=X line tells which ellipse to test against
    side = 0 if p[1] < p[0] else 1
    if inside(p, side):
        cur = side                    # firmly inside a state's ellipse
    # else: ambiguous -> hold previous state
    assigned[k] = cur

# Robust transitions: change in the ellipse-based assigned state
trans_steps = np.where(np.diff(assigned) != 0)[0] + 1
print(f"Total steps: {nsteps}")
print(f"Ellipse-based transition count: {len(trans_steps)}")
for i, step in enumerate(trans_steps):
    frm = assigned[step-1]
    to  = assigned[step]
    print(f"  transition {i+1}: step {step} (t = {step*dt:.2f})  state {frm} -> {to}")

# ----------------------------------------------------------------------
# Crude counter (the naive baseline): label every point purely by which
# side of Y = X it is on, and count every side change as a transition.
# This does NOT use the ellipses, so boundary jitter is counted repeatedly.
# ----------------------------------------------------------------------
crude_state = (Y >= X).astype(int)                 # side of the Y=X line
crude_trans = np.where(np.diff(crude_state) != 0)[0] + 1
print(f"Crude (Y=X side) transition count: {len(crude_trans)}")

# Show that a single crossing event produces a CLUSTER of crude transitions
# within a few steps (jitter), whereas the robust count collapses them to one.
if len(crude_trans) > 1:
    gaps = np.diff(crude_trans)
    clusters = np.sum(gaps > 5) + 1                # crossings separated by >5 steps
    print(f"Crude transitions separated by >5 steps (event clusters): {clusters}")
    print(f"Crude count / robust count ratio: {len(crude_trans)/max(len(trans_steps),1):.2f}")
    # illustrate one cluster explicitly
    example = crude_trans[:min(6, len(crude_trans))]
    print(f"Example early crude transition steps (note tight clustering): {example.tolist()}")

# ----------------------------------------------------------------------
# Phase plane with the two fitted (enlarged) state ellipses.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 7))
ax.plot(X, Y, lw=0.3, color="0.6", alpha=0.6, zorder=1)
colors = {0: "tab:blue", 1: "tab:red"}
for s, cloud in clouds.items():
    ax.scatter(cloud[:, 0], cloud[:, 1], s=2, color=colors[s], alpha=0.15, zorder=2)
for mu, w, h, angle, s in ell_patches:
    e = Ellipse(mu, w, h, angle=angle, fill=False,
                edgecolor=colors[s], lw=2.5, zorder=4,
                label=f"state {s} ellipse (1.5x)")
    ax.add_patch(e)
    ax.plot(*mu, "o", color=colors[s], zorder=5)
lim = [0, max(X.max(), Y.max()) * 1.05]
ax.plot(lim, lim, "k--", lw=1, alpha=0.7, label="split line Y = X")
ax.set_xlim(lim); ax.set_ylim(lim)
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Toggle-switch phase plane with fitted state ellipses (1.5x)")
ax.legend(loc="upper right", fontsize=8)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.2.1_s3.png", dpi=130)

# ----------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms the result because a single genuine "
      "state jump makes the crude side-of-line counter fire many times within "
      "a few steps (the trajectory jitters back and forth across Y=X), so the "
      "crude count vastly exceeds the ellipse-based count, showing the crude "
      "statistic over-counts and the ellipse method's hysteresis is needed.")
