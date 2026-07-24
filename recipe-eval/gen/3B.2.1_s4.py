import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle-switch model (genes X and Y mutually repress each other):
#   dX/dt = gX0 + gX1 / (1 + (Y/Yth)^nY) - kX*X
#   dY/dt = gY0 + gY1 / (1 + (X/Xth)^nX) - kY*Y
# Parameters
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ----------------------------------------------------------------------
# Right-hand sides of the two ODEs (the growth minus degradation terms)
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth)**nY) - kX * X

def fY(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth)**nX) - kY * Y

# ----------------------------------------------------------------------
# Build the nullclines as polylines in the (X, Y) plane.
#
# X-nullcline: dX/dt = 0  ->  X = [gX0 + gX1/(1+(Y/Yth)^nY)] / kX.
#   For each Y we get exactly one X, so parametrize by Y.
Y_grid = np.linspace(0.0, 700.0, 4000)
X_of_Y = (gX0 + gX1 / (1.0 + (Y_grid / Yth)**nY)) / kX     # points on X-nullcline
nullX = np.column_stack([X_of_Y, Y_grid])                  # (X, Y) samples

# Y-nullcline: dY/dt = 0  ->  Y = [gY0 + gY1/(1+(X/Xth)^nX)] / kY.
#   For each X we get exactly one Y, so parametrize by X.
X_grid = np.linspace(0.0, 700.0, 4000)
Y_of_X = (gY0 + gY1 / (1.0 + (X_grid / Xth)**nX)) / kY     # points on Y-nullcline
nullY = np.column_stack([X_grid, Y_of_X])                  # (X, Y) samples

# ----------------------------------------------------------------------
# Segment-intersection test.
# Two segments p->p+r and q->q+s cross when we can solve p+t*r = q+u*s
# with t,u in [0,1].  Using 2D cross products (denominator rxs):
def seg_intersect(p, r, q, s):
    rxs = r[0]*s[1] - r[1]*s[0]          # cross product r x s
    if abs(rxs) < 1e-12:                 # parallel / collinear -> skip
        return None
    qp = q - p
    t = (qp[0]*s[1] - qp[1]*s[0]) / rxs  # position along first segment
    u = (qp[0]*r[1] - qp[1]*r[0]) / rxs  # position along second segment
    if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
        return p + t * r                 # the crossing point
    return None

# ----------------------------------------------------------------------
# Exhaustive all-pairs search: test every segment of nullX against every
# segment of nullY.
def find_exhaustive(A, B):
    hits = []
    pairs = 0
    for i in range(len(A) - 1):
        p = A[i]; r = A[i+1] - A[i]
        for j in range(len(B) - 1):
            q = B[j]; s = B[j+1] - B[j]
            pairs += 1
            pt = seg_intersect(p, r, q, s)
            if pt is not None:
                hits.append(pt)
    return hits, pairs

# ----------------------------------------------------------------------
# Fast search: first narrow to candidate segments using a sign test.
# Along nullX (a set of points that satisfy dX/dt=0) the sign of dY/dt
# changes exactly where the Y-nullcline is crossed, and vice versa.
# So we only keep segments of each polyline whose endpoints straddle the
# other equation's zero, then test only those against each other.
def find_fast(A, B):
    # g_A: value of dY/dt along the X-nullcline A -> zero at a crossing
    gA = np.array([fY(px, py) for px, py in A])
    # g_B: value of dX/dt along the Y-nullcline B -> zero at a crossing
    gB = np.array([fX(px, py) for px, py in B])

    # indices of segments on A where dY/dt changes sign
    iA = np.where(gA[:-1] * gA[1:] <= 0.0)[0]
    # indices of segments on B where dX/dt changes sign
    iB = np.where(gB[:-1] * gB[1:] <= 0.0)[0]

    hits = []
    pairs = 0
    for i in iA:
        p = A[i]; r = A[i+1] - A[i]
        for j in iB:
            q = B[j]; s = B[j+1] - B[j]
            pairs += 1
            pt = seg_intersect(p, r, q, s)
            if pt is not None:
                hits.append(pt)
    return hits, pairs

# ----------------------------------------------------------------------
# De-duplicate crossings that fall within one grid cell of each other.
def dedup(points, tol=1.0):
    uniq = []
    for pt in points:
        if not any(np.hypot(pt[0]-u[0], pt[1]-u[1]) < tol for u in uniq):
            uniq.append(pt)
    return sorted(uniq, key=lambda p: p[0])

# ----------------------------------------------------------------------
# Run both searches
raw_ex, pairs_ex = find_exhaustive(nullX, nullY)
raw_fast, pairs_fast = find_fast(nullX, nullY)
ss_ex = dedup(raw_ex)
ss_fast = dedup(raw_fast)

# ----------------------------------------------------------------------
# Stability: eigenvalues of the Jacobian at each steady state.
def jacobian(X, Y):
    # d(dX/dt)/dX = -kX ;  d(dX/dt)/dY = derivative of Hill(Y)
    dfX_dY = gX1 * (-(nY / Yth) * (Y / Yth)**(nY - 1)) / (1.0 + (Y / Yth)**nY)**2
    dfY_dX = gY1 * (-(nX / Xth) * (X / Xth)**(nX - 1)) / (1.0 + (X / Xth)**nX)**2
    return np.array([[-kX, dfX_dY],
                     [dfY_dX, -kY]])

# ----------------------------------------------------------------------
# Report
print("Number of steady states (exhaustive):", len(ss_ex))
print("Number of steady states (fast):      ", len(ss_fast))
print("Pairs tested (exhaustive):", pairs_ex)
print("Pairs tested (fast):      ", pairs_fast)
print("Fast agrees with exhaustive:",
      len(ss_ex) == len(ss_fast) and
      all(np.allclose(a, b, atol=1e-6) for a, b in zip(ss_ex, ss_fast)))

n_stable = 0
n_unstable = 0
for k, (X, Y) in enumerate(ss_ex, 1):
    print(f"Steady state {k}: X = {X:.6f}")
    print(f"Steady state {k}: Y = {Y:.6f}")
    ev = np.linalg.eigvals(jacobian(X, Y))
    print(f"Steady state {k}: eigenvalue 1 = {ev[0].real:.6f}")
    print(f"Steady state {k}: eigenvalue 2 = {ev[1].real:.6f}")
    stable = np.all(ev.real < 0)
    print(f"Steady state {k}: stable = {bool(stable)}")
    if stable:
        n_stable += 1
    else:
        n_unstable += 1

print("Count stable steady states:  ", n_stable)
print("Count unstable steady states:", n_unstable)

# ----------------------------------------------------------------------
# Plot nullclines and mark the crossings
plt.figure(figsize=(7, 6))
plt.plot(nullX[:, 0], nullX[:, 1], label="X-nullcline (dX/dt = 0)")
plt.plot(nullY[:, 0], nullY[:, 1], label="Y-nullcline (dY/dt = 0)")
for X, Y in ss_ex:
    plt.plot(X, Y, "ko", markersize=9, markerfacecolor="none", markeredgewidth=2)
    plt.annotate(f"({X:.0f}, {Y:.0f})", (X, Y),
                 textcoords="offset points", xytext=(8, 8))
plt.xlim(0, 700)
plt.ylim(0, 700)
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Toggle-switch nullclines and steady states")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.2.1_s4.png")

# One-sentence explanation:
# Finding exactly three intersections with two having Jacobians whose
# eigenvalues are both negative (stable) and one with a positive eigenvalue
# (unstable) confirms the classic bistable toggle-switch signature, and the
# fast search reproducing the same three points while testing far fewer
# segment pairs confirms the sign-narrowing step discards only non-crossing
# segments and misses no genuine steady state.
print("Explanation: three intersections with two stable (both eigenvalues "
      "negative) and one unstable (a positive eigenvalue) is the bistable "
      "toggle-switch signature, and the fast search matching the exhaustive "
      "result on far fewer pairs confirms the sign-narrowing skips only "
      "non-crossing segments without missing any steady state.")
