import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 6D.1 model: noisy genetic toggle switch (two mutually repressing genes).
#   dX = ( b/(1+Y^n) - X ) dt + sigma * dW_x
#   dY = ( b/(1+X^n) - Y ) dt + sigma * dW_y
# Same trajectory as 6D.1: seed 3, b = 20, dt = 0.01, integrated to t = 1000.
# ----------------------------------------------------------------------
np.random.seed(3)
b      = 20.0
n      = 2               # Hill coefficient
sigma  = 6.0             # noise strength (drives the switching)
dt     = 0.01
T      = 1000.0
N      = int(round(T / dt))
sqdt   = np.sqrt(dt)

X = np.empty(N)
Y = np.empty(N)
X[0], Y[0] = 20.0, 0.05  # start in the high-X state

# Euler-Maruyama integration of the SDE
for i in range(N - 1):
    fx = b / (1.0 + Y[i]**n) - X[i]          # drift on X
    fy = b / (1.0 + X[i]**n) - Y[i]          # drift on Y
    X[i+1] = X[i] + fx * dt + sigma * sqdt * np.random.randn()
    Y[i+1] = Y[i] + fy * dt + sigma * sqdt * np.random.randn()

P = np.column_stack([X, Y])   # (N,2) points in the phase plane

print("=== Trajectory (6D.1) ===")
print(f"seed=3  b={b}  n={n}  sigma={sigma}  dt={dt}  T={T}  steps N={N}")

# ----------------------------------------------------------------------
# Split the plane along the line Y = X into the two state clouds.
#   State A = high-X cloud  (X > Y)
#   State B = high-Y cloud  (Y > X)
# ----------------------------------------------------------------------
maskA = X > Y
maskB = Y > X
ENLARGE = 1.5  # factor by which each ellipse's axes are enlarged

def fit_state(points):
    """Fit an ellipse to a cloud via covariance eigen-decomposition."""
    mu   = points.mean(axis=0)                 # ellipse center
    C    = np.cov(points.T)                     # covariance matrix
    lam, V = np.linalg.eigh(C)                  # eigenvalues/eigenvectors
    # half-axis lengths at 1-sigma, then enlarged by ENLARGE
    half_axes = ENLARGE * np.sqrt(lam)
    Cinv = np.linalg.inv(C)
    return mu, C, lam, V, half_axes, Cinv

muA, CA, lamA, VA, axA, CinvA = fit_state(P[maskA])
muB, CB, lamB, VB, axB, CinvB = fit_state(P[maskB])

def inside(points, mu, Cinv, thresh=ENLARGE):
    """Point is inside the enlarged ellipse when Mahalanobis distance <= thresh."""
    d = points - mu
    m2 = np.einsum('ij,jk,ik->i', d, Cinv, d)  # squared Mahalanobis distance
    return np.sqrt(m2) <= thresh

inA = inside(P, muA, CinvA)   # inside state-A ellipse
inB = inside(P, muB, CinvB)   # inside state-B ellipse

print("\n=== Fitted state ellipses (axes enlarged x1.5) ===")
print(f"State A center = ({muA[0]:.3f}, {muA[1]:.3f})   half-axes = "
      f"({axA[0]:.3f}, {axA[1]:.3f})   points in cloud = {maskA.sum()}")
print(f"State B center = ({muB[0]:.3f}, {muB[1]:.3f})   half-axes = "
      f"({axB[0]:.3f}, {axB[1]:.3f})   points in cloud = {maskB.sum()}")

# ----------------------------------------------------------------------
# Robust (ellipse) transition counting.
# Assign each point to a state ONLY when it lies inside exactly one ellipse;
# in the "dead zone" (inside both, or inside neither) carry the previous
# label forward.  A transition is a change of assigned label.
# ----------------------------------------------------------------------
labels = np.empty(N, dtype='<U1')
cur = 'A' if X[0] > Y[0] else 'B'   # initial label from the starting side
for i in range(N):
    if inA[i] and not inB[i]:
        cur = 'A'                    # clearly in state A
    elif inB[i] and not inA[i]:
        cur = 'B'                    # clearly in state B
    # else: ambiguous -> keep previous label (this suppresses jitter)
    labels[i] = cur

# detect transitions from the change in assigned state
ell_steps = np.where(labels[1:] != labels[:-1])[0] + 1

print("\n=== Robust (ellipse) transitions ===")
print(f"Total detected transitions = {len(ell_steps)}")
for k, s in enumerate(ell_steps, 1):
    print(f"  transition {k:3d}:  step = {s:6d}   time = {s*dt:8.2f}   "
          f"{labels[s-1]} -> {labels[s]}")

# ----------------------------------------------------------------------
# Crude counter (the naive baseline).
# It simply counts every sign flip of (X - Y) about the Y = X line.
# ----------------------------------------------------------------------
sgn = np.sign(X - Y)
crude_steps = np.where(np.diff(sgn) != 0)[0] + 1
print("\n=== Crude counter (sign flips of X - Y) ===")
print(f"Total crude transitions = {len(crude_steps)}")

# Show the over-counting: group crude flips that occur within a few steps of
# each other into clusters; a single physical jump produces one such cluster.
WIN = 5
clusters = []
if len(crude_steps):
    start = prev = crude_steps[0]
    members = [prev]
    for s in crude_steps[1:]:
        if s - prev <= WIN:
            members.append(s)          # same boundary-jitter event
        else:
            clusters.append(members)
            members = [s]
        prev = s
    clusters.append(members)

sizes = [len(c) for c in clusters]
biggest = clusters[int(np.argmax(sizes))]
print(f"Crude flips group into {len(clusters)} clusters (gap<= {WIN} steps).")
print(f"Number of clusters with >1 flip (jitter events) = "
      f"{sum(1 for c in clusters if len(c) > 1)}")
print(f"Largest boundary-jitter cluster: {len(biggest)} crude transitions "
      f"within steps {biggest[0]}..{biggest[-1]} "
      f"(span {biggest[-1]-biggest[0]} steps): {biggest}")

# ----------------------------------------------------------------------
# Phase plane with the two fitted state ellipses.
# ----------------------------------------------------------------------
def ellipse_pts(mu, V, half_axes, m=300):
    t = np.linspace(0, 2*np.pi, m)
    unit = np.stack([np.cos(t), np.sin(t)])         # unit circle
    return mu[:, None] + V @ (np.diag(half_axes) @ unit)

eA = ellipse_pts(muA, VA, axA)
eB = ellipse_pts(muB, VB, axB)

fig, ax = plt.subplots(figsize=(7, 7))
sub = slice(None, None, 20)   # subsample for a lighter scatter
ax.scatter(X[sub][maskA[sub]], Y[sub][maskA[sub]], s=3, c='tab:blue',
           alpha=0.3, label='X > Y (state A)')
ax.scatter(X[sub][maskB[sub]], Y[sub][maskB[sub]], s=3, c='tab:red',
           alpha=0.3, label='Y > X (state B)')
ax.plot(eA[0], eA[1], 'b-', lw=2, label='State A ellipse (x1.5)')
ax.plot(eB[0], eB[1], 'r-', lw=2, label='State B ellipse (x1.5)')
lim = [min(X.min(), Y.min()), max(X.max(), Y.max())]
ax.plot(lim, lim, 'k--', lw=1, label='Y = X split')
ax.plot(*muA, 'b*', ms=14)
ax.plot(*muB, 'r*', ms=14)
ax.set_xlabel('X'); ax.set_ylabel('Y')
ax.set_title('Toggle-switch phase plane with fitted state ellipses')
ax.legend(loc='upper right', fontsize=8)
ax.set_aspect('equal', 'box')
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.2.1_s1.png")

# ----------------------------------------------------------------------
print("\n=== Summary ===")
print(f"Robust (ellipse) transition count = {len(ell_steps)}")
print(f"Crude (sign-flip)  transition count = {len(crude_steps)}")
print("Explanation: The crude counter over-counts because a single physical "
      "switch lingers on the Y=X boundary and flips the sign of X-Y many times "
      "within a few steps, producing a tight cluster of crude transitions, "
      "whereas the dead zone between the two enlarged ellipses forces the point "
      "to actually reach the opposite cloud before the assigned state changes, "
      "collapsing each such cluster into one robust transition and confirming "
      "the ellipse count is the reliable statistic.")
