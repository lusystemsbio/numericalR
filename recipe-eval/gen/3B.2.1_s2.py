import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Right-hand sides of the ODEs
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def fY(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ----------------------------------------------------------------------
# Build the two nullclines as polylines (dense point lists).
# X-nullcline (dX/dt = 0):  X is an explicit function of Y.
# Y-nullcline (dY/dt = 0):  Y is an explicit function of X.
# ----------------------------------------------------------------------
Ygrid = np.linspace(0.0, 400.0, 4000)                       # parameter for X-nullcline
Xnull_X = (gX0 + gX1 / (1.0 + (Ygrid / Yth) ** nY)) / kX     # X value on X-nullcline
Xnull = np.column_stack([Xnull_X, Ygrid])                   # points (X, Y)

Xgrid = np.linspace(0.0, 600.0, 4000)                       # parameter for Y-nullcline
Ynull_Y = (gY0 + gY1 / (1.0 + (Xgrid / Xth) ** nX)) / kY     # Y value on Y-nullcline
Ynull = np.column_stack([Xgrid, Ynull_Y])                   # points (X, Y)

# ----------------------------------------------------------------------
# Segment-intersection helpers (explicit, no library routine).
# ----------------------------------------------------------------------
def orient(a, b, c):
    # Signed area sign of triangle (a, b, c): >0 CCW, <0 CW, 0 collinear.
    return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

def segments_cross(p1, p2, p3, p4):
    # Proper segment-segment intersection test via orientation signs.
    d1 = orient(p3, p4, p1)
    d2 = orient(p3, p4, p2)
    d3 = orient(p1, p2, p3)
    d4 = orient(p1, p2, p4)
    return (((d1 > 0) != (d2 > 0)) and ((d3 > 0) != (d4 > 0)))

def intersection_point(p1, p2, p3, p4):
    # Line-line intersection of the two segments (they are known to cross).
    x1, y1 = p1; x2, y2 = p2; x3, y3 = p3; x4, y4 = p4
    den = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    px = ((x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)) / den
    py = ((x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)) / den
    return (px, py)

# ----------------------------------------------------------------------
# Exhaustive all-pairs search: test every segment of curve A against
# every segment of curve B.
# ----------------------------------------------------------------------
def find_intersections_exhaustive(A, B):
    hits = []
    tested = 0
    for i in range(len(A) - 1):
        for j in range(len(B) - 1):
            tested += 1                                  # count every pair examined
            if segments_cross(A[i], A[i + 1], B[j], B[j + 1]):
                hits.append(intersection_point(A[i], A[i + 1], B[j], B[j + 1]))
    return hits, tested

# ----------------------------------------------------------------------
# Fast search: first narrow to segments whose bounding boxes overlap.
# A pair can only cross if their x-intervals AND y-intervals overlap;
# this is a cheap sign/interval test that rejects almost all pairs
# before the (more expensive) orientation-based crossing test runs.
# ----------------------------------------------------------------------
def find_intersections_fast(A, B):
    hits = []
    tested = 0
    # Precompute per-segment axis-aligned bounding boxes for curve B.
    Bmin_x = np.minimum(B[:-1, 0], B[1:, 0]); Bmax_x = np.maximum(B[:-1, 0], B[1:, 0])
    Bmin_y = np.minimum(B[:-1, 1], B[1:, 1]); Bmax_y = np.maximum(B[:-1, 1], B[1:, 1])
    for i in range(len(A) - 1):
        ax0, ax1 = A[i, 0], A[i + 1, 0]
        ay0, ay1 = A[i, 1], A[i + 1, 1]
        axmin, axmax = min(ax0, ax1), max(ax0, ax1)
        aymin, aymax = min(ay0, ay1), max(ay0, ay1)
        # Vectorized overlap test picks only the candidate segments of B.
        overlap = (Bmax_x >= axmin) & (Bmin_x <= axmax) & \
                  (Bmax_y >= aymin) & (Bmin_y <= aymax)
        for j in np.nonzero(overlap)[0]:
            tested += 1                                  # only overlapping pairs tested
            if segments_cross(A[i], A[i + 1], B[j], B[j + 1]):
                hits.append(intersection_point(A[i], A[i + 1], B[j], B[j + 1]))
    return hits, tested

# Deduplicate nearly-identical crossings.
def dedup(points, tol=1e-3):
    out = []
    for p in points:
        if not any(abs(p[0] - q[0]) < tol and abs(p[1] - q[1]) < tol for q in out):
            out.append(p)
    return sorted(out)

# ----------------------------------------------------------------------
# Run both searches.
# ----------------------------------------------------------------------
hits_ex, tested_ex = find_intersections_exhaustive(Xnull, Ynull)
hits_fa, tested_fa = find_intersections_fast(Xnull, Ynull)

ss_ex = dedup(hits_ex)
ss_fa = dedup(hits_fa)

print("Pairs tested (exhaustive all-pairs):", tested_ex)
print("Pairs tested (fast, sign/box-narrowed):", tested_fa)
print("Number of steady states (exhaustive):", len(ss_ex))
print("Number of steady states (fast):", len(ss_fa))

agree = (len(ss_ex) == len(ss_fa)) and all(
    abs(a[0] - b[0]) < 1e-2 and abs(a[1] - b[1]) < 1e-2 for a, b in zip(ss_ex, ss_fa))
print("Fast search agrees with exhaustive search:", agree)

# ----------------------------------------------------------------------
# Stability of each steady state via the Jacobian eigenvalues.
# ----------------------------------------------------------------------
def jacobian(X, Y):
    dfX_dY = -gX1 * nY * (Y / Yth) ** (nY - 1) / Yth / (1.0 + (Y / Yth) ** nY) ** 2
    dfY_dX = -gY1 * nX * (X / Xth) ** (nX - 1) / Xth / (1.0 + (X / Xth) ** nX) ** 2
    return np.array([[-kX, dfX_dY], [dfY_dX, -kY]])

n_stable, n_unstable = 0, 0
labels = []
for k, (X, Y) in enumerate(ss_ex, 1):
    eig = np.linalg.eigvals(jacobian(X, Y))
    stable = np.all(eig.real < 0)
    kind = "stable" if stable else "unstable"
    if stable:
        n_stable += 1
    else:
        n_unstable += 1
    labels.append(kind)
    print(f"Steady state {k}: X = {X:.4f}, Y = {Y:.4f} | "
          f"eigenvalues = {eig[0]:.5f}, {eig[1]:.5f} | {kind}")

print("Count stable steady states:", n_stable)
print("Count unstable steady states:", n_unstable)

# ----------------------------------------------------------------------
# Plot the nullclines and mark the crossings (steady states).
# ----------------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(Xnull[:, 0], Xnull[:, 1], label="X-nullcline (dX/dt = 0)", color="tab:blue")
plt.plot(Ynull[:, 0], Ynull[:, 1], label="Y-nullcline (dY/dt = 0)", color="tab:red")
for (X, Y), kind in zip(ss_ex, labels):
    marker = "o" if kind == "stable" else "s"
    plt.scatter([X], [Y], s=120, zorder=5, edgecolor="k",
                color=("green" if kind == "stable" else "orange"),
                marker=marker, label=f"{kind} ({X:.1f}, {Y:.1f})")
plt.xlabel("X"); plt.ylabel("Y")
plt.title("Toggle-switch nullclines and steady states")
plt.legend(fontsize=8); plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.2.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("Check rationale: finding exactly three nullcline crossings whose Jacobian "
      "eigenvalues give two stable and one unstable state is the signature of a "
      "bistable toggle switch, and the fast box-narrowed search reproducing the "
      "exhaustive result while testing far fewer pairs confirms the narrowing "
      "discarded only non-intersecting pairs.")
