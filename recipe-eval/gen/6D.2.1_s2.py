import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 6D.1 model: symmetric two-gene toggle-switch SDE (Euler-Maruyama).
# dx = (b/(1+y^n) - x) dt + sigma dW,  dy = (b/(1+x^n) - y) dt + sigma dW
# Same trajectory settings as 6D.1: seed 3, b = 20, dt = 0.01, T = 1000.
# ---------------------------------------------------------------
np.random.seed(3)
b     = 20.0     # production strength
n     = 4        # Hill coefficient
sigma = 6.0      # noise amplitude (drives the switching)
dt    = 0.01
T     = 1000.0
N     = int(T / dt)

x = np.empty(N + 1)
y = np.empty(N + 1)
x[0], y[0] = b, 0.0          # start in the (high-x, low-y) state
sqdt = np.sqrt(dt)
for i in range(N):
    fx = b / (1.0 + y[i]**n) - x[i]
    fy = b / (1.0 + x[i]**n) - y[i]
    x[i+1] = x[i] + fx*dt + sigma*sqdt*np.random.randn()
    y[i+1] = y[i] + fy*dt + sigma*sqdt*np.random.randn()

pts = np.column_stack([x, y])   # (N+1, 2) trajectory in the phase plane

# ---------------------------------------------------------------
# Step 1: split the plane along the diagonal Y = X to get the two
# raw point clouds, one per state.
#   state A = high-x / low-y  -> points with x >  y
#   state B = low-x  / high-y -> points with x <= y
# ---------------------------------------------------------------
maskA = pts[:, 0] >  pts[:, 1]
maskB = pts[:, 0] <= pts[:, 1]
cloudA = pts[maskA]
cloudB = pts[maskB]

# ---------------------------------------------------------------
# Step 2: fit an ellipse to each cloud from the covariance
# eigen-decomposition.  The eigenvectors give the ellipse axes'
# directions and the eigenvalues their (variance) lengths.
# Semi-axes = k * sqrt(eigenvalue); here k = 2 (a 2-sigma ellipse),
# then every axis is enlarged by the requested factor 1.5.
# ---------------------------------------------------------------
K_STD   = 2.0
ENLARGE = 1.5

def fit_ellipse(cloud):
    mean = cloud.mean(axis=0)
    cov  = np.cov(cloud.T)                 # 2x2 covariance
    evals, evecs = np.linalg.eigh(cov)     # eigen-decomposition (axes)
    semi = ENLARGE * K_STD * np.sqrt(evals)   # enlarged semi-axis lengths
    return mean, evecs, semi

meanA, evecA, semiA = fit_ellipse(cloudA)
meanB, evecB, semiB = fit_ellipse(cloudB)

# ---------------------------------------------------------------
# Step 3: a point is "inside" an ellipse when, expressed in the
# ellipse's own axis frame, (u/a)^2 + (v/a')^2 <= 1.
# ---------------------------------------------------------------
def inside(p, mean, evec, semi):
    d = p - mean                 # offset from ellipse centre
    coords = evec.T @ d.T        # project onto the principal axes (columns of evec)
    return ((coords[0]/semi[0])**2 + (coords[1]/semi[1])**2) <= 1.0

inA = inside(pts, meanA, evecA, semiA)
inB = inside(pts, meanB, evecB, semiB)

# ---------------------------------------------------------------
# Step 4: assign each step to a state (0 = A, 1 = B).  When a point
# is inside neither ellipse (in the gap between states) we carry the
# previously assigned state forward, so a transition is only recorded
# once the trajectory has actually entered the opposite ellipse.
# ---------------------------------------------------------------
assigned = np.full(N + 1, -1, dtype=int)
cur = -1
for i in range(N + 1):
    if inA[i]:
        cur = 0
    elif inB[i]:
        cur = 1
    assigned[i] = cur           # -1 only until the first ellipse is entered

# ---------------------------------------------------------------
# Step 5: a transition is a change in the assigned state between two
# consecutive valid (non -1) assignments.
# ---------------------------------------------------------------
robust_steps = []
for i in range(1, N + 1):
    if assigned[i] != -1 and assigned[i-1] != -1 and assigned[i] != assigned[i-1]:
        robust_steps.append(i)

# ---------------------------------------------------------------
# Crude counter (the over-counting baseline): count every crossing
# of the Y = X boundary, i.e. every sign change of (x - y).
# ---------------------------------------------------------------
diff = pts[:, 0] - pts[:, 1]
sign = np.sign(diff)
crude_steps = [i for i in range(1, N + 1) if sign[i] != sign[i-1] and sign[i] != 0]

# ---------------------------------------------------------------
# Report the robust (ellipse-based) transitions.
# ---------------------------------------------------------------
print(f"Ellipse-based (robust) transition count: {len(robust_steps)}")
for k, s in enumerate(robust_steps, 1):
    frm = 'A' if assigned[s-1] == 0 else 'B'
    to  = 'A' if assigned[s]   == 0 else 'B'
    print(f"  transition {k}: step {s} (t = {s*dt:.2f}), {frm} -> {to}")

# ---------------------------------------------------------------
# Check: show the crude counter over-counts.  A single physical
# boundary-jitter event produces a burst of Y=X crossings within a
# few steps; we flag crude transitions whose gap to the previous one
# is <= 5 steps as members of such jitter clusters.
# ---------------------------------------------------------------
print(f"\nCrude (Y=X sign-change) transition count: {len(crude_steps)}")
crude_arr = np.array(crude_steps)
gaps = np.diff(crude_arr)
jitter_members = int(np.sum(gaps <= 5))
print(f"Crude transitions separated by <= 5 steps (jitter-cluster members): {jitter_members}")
print(f"Over-count factor crude/robust: {len(crude_steps) / max(len(robust_steps),1):.1f}x")

# Illustrate one cluster explicitly: the tightest burst of crude crossings.
if len(gaps) > 0:
    j = int(np.argmin(gaps))
    lo = j
    while lo > 0 and gaps[lo-1] <= 5:
        lo -= 1
    hi = j
    while hi < len(gaps) and gaps[hi] <= 5:
        hi += 1
    cluster = crude_arr[lo:hi+1]
    print(f"Example jitter cluster (steps): {cluster.tolist()}  "
          f"-> {len(cluster)} crude crossings inside a {cluster[-1]-cluster[0]}-step window")

# ---------------------------------------------------------------
# Plot: phase plane with the two fitted (enlarged) state ellipses.
# ---------------------------------------------------------------
def ellipse_curve(mean, evec, semi):
    th = np.linspace(0, 2*np.pi, 200)
    circ = np.vstack([semi[0]*np.cos(th), semi[1]*np.sin(th)])  # in axis frame
    return (mean[:, None] + evec @ circ)                        # back to plane

fig, ax = plt.subplots(figsize=(7, 7))
ax.plot(pts[::20, 0], pts[::20, 1], '.', ms=1, color='0.7', alpha=0.5, label='trajectory')
lim = [min(pts.min(), -2), pts.max() + 2]
ax.plot(lim, lim, 'k--', lw=1, label='Y = X split')

eA = ellipse_curve(meanA, evecA, semiA)
eB = ellipse_curve(meanB, evecB, semiB)
ax.plot(eA[0], eA[1], 'r-', lw=2, label='state A ellipse (x > y)')
ax.plot(eB[0], eB[1], 'b-', lw=2, label='state B ellipse (x < y)')
ax.plot(*meanA, 'r+', ms=12, mew=2)
ax.plot(*meanB, 'b+', ms=12, mew=2)

ax.set_xlabel('X'); ax.set_ylabel('Y')
ax.set_title('Toggle switch: fitted state ellipses (axes x1.5)')
ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
ax.legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.2.1_s2.png", dpi=120)

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ---------------------------------------------------------------
print("\nWhy the check confirms the result: because the crude counter fires many "
      "times within a few-step window around a single Y=X crossing (a jitter cluster) "
      "while the ellipse method requires the trajectory to actually enter the opposite "
      "state's cloud, the much larger crude count exposes the over-counting and confirms "
      "that the ellipse-based statistic reports genuine state switches rather than boundary noise.")
