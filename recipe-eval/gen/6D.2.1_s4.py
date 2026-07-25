import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------
# 1. Reproduce the 6D.1 noisy toggle-switch SDE and trajectory
#    (genetic toggle switch: X and Y mutually repress each other)
#    Euler-Maruyama, seed 3, b = 20, dt = 0.01, up to t = 1000.
# ------------------------------------------------------------------
seed = 3
b = 20.0
n = 2          # Hill coefficient
sigma = 5.0    # noise amplitude
dt = 0.01
T = 1000.0
nsteps = int(T / dt)

rng = np.random.default_rng(seed)

def drift(x, y):
    # mutual repression: high Y suppresses X production and vice versa
    dx = b / (1.0 + y**n) - x
    dy = b / (1.0 + x**n) - y
    return dx, dy

X = np.empty(nsteps + 1)
Y = np.empty(nsteps + 1)
X[0], Y[0] = 20.0, 0.0          # start in the X-dominant state
sqdt = np.sqrt(dt)
for i in range(nsteps):
    dx, dy = drift(X[i], Y[i])
    X[i + 1] = X[i] + dx * dt + sigma * sqdt * rng.standard_normal()
    Y[i + 1] = Y[i] + dy * dt + sigma * sqdt * rng.standard_normal()

pts = np.column_stack([X, Y])   # shape (N, 2)

# ------------------------------------------------------------------
# 2. Split the plane along Y = X to get the two state clouds
#    State 0 : X-dominant  (X > Y)
#    State 1 : Y-dominant  (Y > X)
# ------------------------------------------------------------------
maskA = pts[:, 0] >= pts[:, 1]   # cloud for state 0
maskB = ~maskA                   # cloud for state 1
cloudA = pts[maskA]
cloudB = pts[maskB]

# ------------------------------------------------------------------
# 3. Fit an ellipse to each cloud from the covariance eigen-decomposition.
#    mean = centre, eigenvectors = ellipse axis directions,
#    sqrt(eigenvalue) = axis half-lengths (scaled to a 2-sigma ellipse),
#    then enlarge every axis by the factor 1.5.
# ------------------------------------------------------------------
k = 2.0            # base ellipse at 2 standard deviations
enlarge = 1.5      # enlarge the axes by this factor

def fit_ellipse(cloud):
    mean = cloud.mean(axis=0)
    cov = np.cov(cloud, rowvar=False)               # 2x2 covariance
    eigval, eigvec = np.linalg.eigh(cov)            # eigen-decomposition
    axes = enlarge * k * np.sqrt(eigval)            # enlarged half-axis lengths
    return mean, eigvec, axes

meanA, vecA, axesA = fit_ellipse(cloudA)
meanB, vecB, axesB = fit_ellipse(cloudB)

def inside(p, mean, eigvec, axes):
    # project the centred point onto the ellipse axes and test (u/a)^2+(v/b)^2<=1
    d = p - mean
    proj = eigvec.T @ d                             # coordinates in axis frame
    return np.sum((proj / axes) ** 2) <= 1.0

def maha(p, mean, eigvec, axes):
    # squared normalised distance, used only to break ties / assign leftovers
    d = p - mean
    proj = eigvec.T @ d
    return np.sum((proj / axes) ** 2)

# ------------------------------------------------------------------
# 4. Assign each point to a state (inside which ellipse it lies).
#    inside A only -> 0 ; inside B only -> 1 ;
#    inside both or neither -> whichever ellipse is nearer (Mahalanobis).
# ------------------------------------------------------------------
assigned = np.empty(nsteps + 1, dtype=int)
for i, p in enumerate(pts):
    inA = inside(p, meanA, vecA, axesA)
    inB = inside(p, meanB, vecB, axesB)
    if inA and not inB:
        assigned[i] = 0
    elif inB and not inA:
        assigned[i] = 1
    else:
        # ambiguous: pick the closer ellipse
        assigned[i] = 0 if maha(p, meanA, vecA, axesA) <= maha(p, meanB, vecB, axesB) else 1

# ------------------------------------------------------------------
# 5. Detect a transition wherever the assigned state changes step-to-step
#    (this is the CRUDE per-step-change counter).
# ------------------------------------------------------------------
trans_steps = []
trans_dir = []
for i in range(1, nsteps + 1):
    if assigned[i] != assigned[i - 1]:
        trans_steps.append(i)
        trans_dir.append((assigned[i - 1], assigned[i]))
trans_steps = np.array(trans_steps, dtype=int)

crude_count = len(trans_steps)

names = {0: "X-dominant(state0)", 1: "Y-dominant(state1)"}
print("=== Detected transitions (crude per-step-change counter) ===")
for j, s in enumerate(trans_steps):
    a, bstate = trans_dir[j]
    print(f"Transition {j+1}: step = {s}: {names[a]} -> {names[bstate]}")
print(f"Crude transition count = {crude_count}")

# ------------------------------------------------------------------
# 6. Cross-check: group transitions that fall within a few steps of each
#    other into clusters. A single true switch should be one isolated
#    event; clusters reveal boundary jitter that the crude counter inflates.
# ------------------------------------------------------------------
gap = 20   # "a few steps" tolerance for grouping (=0.2 time units)
clusters = []
if crude_count > 0:
    current = [trans_steps[0]]
    for s in trans_steps[1:]:
        if s - current[-1] <= gap:
            current.append(s)          # same jitter cluster
        else:
            clusters.append(current)   # start a new cluster
            current = [s]
    clusters.append(current)

robust_count = len(clusters)

print("\n=== Cross-check: clustering of transitions ===")
for c, cl in enumerate(clusters):
    print(f"Cluster {c+1}: {len(cl)} crude transitions between steps {cl[0]} and {cl[-1]}")
multi = [cl for cl in clusters if len(cl) > 1]
print(f"Number of clusters (robust switch estimate) = {robust_count}")
print(f"Crude count = {crude_count}")
print(f"Clusters containing more than one crude transition (jitter bursts) = {len(multi)}")
print(f"Over-count factor (crude / robust) = {crude_count / robust_count if robust_count else float('nan')}")
print("Crude counter over-counts:", crude_count > robust_count)

# ------------------------------------------------------------------
# 7. Phase plane with the two fitted (enlarged) state ellipses.
# ------------------------------------------------------------------
theta = np.linspace(0, 2 * np.pi, 200)

def ellipse_curve(mean, eigvec, axes):
    circ = np.column_stack([axes[0] * np.cos(theta), axes[1] * np.sin(theta)])
    return (eigvec @ circ.T).T + mean

eA = ellipse_curve(meanA, vecA, axesA)
eB = ellipse_curve(meanB, vecB, axesB)

fig, ax = plt.subplots(figsize=(7, 7))
ax.scatter(cloudA[:, 0], cloudA[:, 1], s=1, color="tab:blue", alpha=0.15, label="X-dominant cloud")
ax.scatter(cloudB[:, 0], cloudB[:, 1], s=1, color="tab:red", alpha=0.15, label="Y-dominant cloud")
ax.plot(eA[:, 0], eA[:, 1], color="navy", lw=2, label="state-0 ellipse (x1.5)")
ax.plot(eB[:, 0], eB[:, 1], color="darkred", lw=2, label="state-1 ellipse (x1.5)")
lim = [min(X.min(), Y.min()) - 2, max(X.max(), Y.max()) + 2]
ax.plot(lim, lim, "k--", lw=1, label="Y = X split")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Toggle-switch phase plane with fitted state ellipses")
ax.set_aspect("equal")
ax.legend(loc="upper right", markerscale=6, fontsize=8)
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.2.1_s4.png")

# ------------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# Because genuine state switches are physically rare and widely separated
# in time, the many transitions occurring within a few steps of each other
# can only be one trajectory wobbling back and forth across the Y=X boundary,
# which proves the crude per-step-change counter inflates a single physical
# switch into a whole cluster and therefore over-counts.
# ------------------------------------------------------------------
print("\nExplanation: genuine switches are rare and well separated, so transitions"
      " bunched within a few steps must come from one boundary-crossing wobbling"
      " in and out of an ellipse, confirming the crude counter over-counts.")
