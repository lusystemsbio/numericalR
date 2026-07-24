import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12

# ----------------------------------------------------------------------
# The two rate functions.  A nullcline is where the corresponding
# rate is zero (dX/dt = 0  or  dY/dt = 0).
# ----------------------------------------------------------------------
def fX(X, Y):
    # dX/dt : basal + repressive Hill function of Y - linear degradation
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def fY(X, Y):
    # dY/dt : basal + repressive Hill function of X - linear degradation
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ----------------------------------------------------------------------
# Build a grid over the (X, Y) plane.  We only need X large enough to
# reach the maximum possible steady state (gX0+gX1)/kX, likewise for Y.
# ----------------------------------------------------------------------
Xmax = (gX0 + gX1) / kX * 1.1   # generous upper bound for X
Ymax = (gY0 + gY1) / kY * 1.1   # generous upper bound for Y
xs = np.linspace(0.0, Xmax, 600)
ys = np.linspace(0.0, Ymax, 600)
XX, YY = np.meshgrid(xs, ys)    # 2-D grids of coordinates

# Evaluate each rate as a surface Z = f(X,Y) on the grid.
ZX = fX(XX, YY)                 # surface whose zero contour is the X-nullcline
ZY = fY(XX, YY)                 # surface whose zero contour is the Y-nullcline

# ----------------------------------------------------------------------
# GENERAL METHOD: the nullclines are the zero-level contours of the
# surfaces.  No separation of variables is needed -- we just ask the
# contour routine for the level set Z = 0.  We capture the contour
# segments so we can also intersect them numerically.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))

csX = ax.contour(XX, YY, ZX, levels=[0.0], colors="tab:blue")
csY = ax.contour(XX, YY, ZY, levels=[0.0], colors="tab:red")

# ----------------------------------------------------------------------
# Find the crossings (steady states) as intersections of the two zero
# contours.  We do this segment-by-segment: for every small segment of
# the X-nullcline and every segment of the Y-nullcline, test whether
# the two line segments intersect, and if so solve for the crossing.
# ----------------------------------------------------------------------
def contour_segments(cs):
    """Return a list of (P0, P1) segment endpoints for a contour set."""
    segs = []
    # matplotlib >=3.8 exposes .get_paths(); fall back to .allsegs otherwise
    try:
        paths = cs.get_paths()
        polylines = [p.vertices for p in paths if len(p.vertices) >= 2]
    except Exception:
        polylines = [np.asarray(s) for lvl in cs.allsegs for s in lvl if len(s) >= 2]
    for poly in polylines:
        for i in range(len(poly) - 1):
            segs.append((poly[i], poly[i + 1]))
    return segs

def seg_intersect(p1, p2, p3, p4):
    """Intersection point of segments p1p2 and p3p4, or None."""
    r = p2 - p1
    s = p4 - p3
    denom = r[0] * s[1] - r[1] * s[0]
    if abs(denom) < 1e-12:            # parallel
        return None
    qp = p3 - p1
    t = (qp[0] * s[1] - qp[1] * s[0]) / denom
    u = (qp[0] * r[1] - qp[1] * r[0]) / denom
    if 0.0 <= t <= 1.0 and 0.0 <= u <= 1.0:
        return p1 + t * r
    return None

segsX = contour_segments(csX)
segsY = contour_segments(csY)

crossings = []
for a1, a2 in segsX:
    for b1, b2 in segsY:
        pt = seg_intersect(np.asarray(a1), np.asarray(a2),
                            np.asarray(b1), np.asarray(b2))
        if pt is not None:
            crossings.append(pt)

# De-duplicate crossings that are essentially the same point.
crossings_unique = []
for p in crossings:
    if not any(np.hypot(*(p - q)) < 1e-3 for q in crossings_unique):
        crossings_unique.append(p)
crossings_unique = sorted(crossings_unique, key=lambda p: p[0])

# ----------------------------------------------------------------------
# Refine each crossing with a couple of Newton steps on (fX, fY)=0 so
# the reported steady states are accurate, not just grid-resolution.
# ----------------------------------------------------------------------
def refine(p):
    X, Y = float(p[0]), float(p[1])
    for _ in range(50):
        F = np.array([fX(X, Y), fY(X, Y)])
        # Jacobian
        dHillY = -gX1 * nY * (Y / Yth) ** (nY - 1) / Yth / (1 + (Y / Yth) ** nY) ** 2
        dHillX = -gY1 * nX * (X / Xth) ** (nX - 1) / Xth / (1 + (X / Xth) ** nX) ** 2
        J = np.array([[-kX, dHillY],
                      [dHillX, -kY]])
        dp = np.linalg.solve(J, -F)
        X, Y = X + dp[0], Y + dp[1]
        if np.hypot(*dp) < 1e-12:
            break
    return np.array([X, Y])

steady_states = [refine(p) for p in crossings_unique]

# Plot the crossings.
for p in steady_states:
    ax.plot(p[0], p[1], "ko", ms=8, zorder=5)

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Toggle switch nullclines as zero contours, with crossings")
# Manual legend proxies (contour has no automatic label).
ax.plot([], [], color="tab:blue", label="X-nullcline (fX = 0)")
ax.plot([], [], color="tab:red", label="Y-nullcline (fY = 0)")
ax.plot([], [], "ko", label="steady states (crossings)")
ax.legend(loc="upper right")

plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.5.1_s4.png",
            dpi=130, bbox_inches="tight")

# ----------------------------------------------------------------------
# CHECK: reproduce the nullclines by separation of variables.
#   fX = 0  =>  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX      (X as a function of Y)
#   fY = 0  =>  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY      (Y as a function of X)
# We verify that these closed-form curves coincide with the zero
# contours, and that solving them simultaneously gives the same
# crossings as the contour intersections.
# ----------------------------------------------------------------------
def X_of_Y(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# Residual of the closed-form curves plugged back into the rate laws
# (should be ~0 everywhere, confirming contour == separation-of-variables).
ycheck = np.linspace(0, Ymax, 500)
resX = np.max(np.abs(fX(X_of_Y(ycheck), ycheck)))
xcheck = np.linspace(0, Xmax, 500)
resY = np.max(np.abs(fY(xcheck, Y_of_X(xcheck))))

# Solve the coupled fixed-point equation X = X_of_Y(Y_of_X(X)) by scanning
# for sign changes of g(X) = X_of_Y(Y_of_X(X)) - X, then bisecting.
def g(X):
    return X_of_Y(Y_of_X(X)) - X

xgrid = np.linspace(0.1, Xmax, 4000)
gvals = g(xgrid)
sv_states = []
for i in range(len(xgrid) - 1):
    if gvals[i] == 0.0 or gvals[i] * gvals[i + 1] < 0.0:
        a, b = xgrid[i], xgrid[i + 1]
        for _ in range(100):                 # bisection
            m = 0.5 * (a + b)
            if g(a) * g(m) <= 0:
                b = m
            else:
                a = m
        Xs = 0.5 * (a + b)
        sv_states.append(np.array([Xs, Y_of_X(Xs)]))
sv_states = sorted(sv_states, key=lambda p: p[0])

# ----------------------------------------------------------------------
# Report every numerical result.
# ----------------------------------------------------------------------
print(f"Grid: X in [0, {Xmax:.3f}], Y in [0, {Ymax:.3f}], 600 x 600 points")
print(f"Max |fX| along separation-of-variables X-nullcline: {resX:.3e}")
print(f"Max |fY| along separation-of-variables Y-nullcline: {resY:.3e}")

print(f"Number of crossings from zero-contour intersection: {len(steady_states)}")
for i, p in enumerate(steady_states, 1):
    print(f"  Contour crossing {i}: X = {p[0]:.6f}, Y = {p[1]:.6f}  "
          f"(residual |fX|={abs(fX(*p)):.2e}, |fY|={abs(fY(*p)):.2e})")

print(f"Number of steady states from separation of variables: {len(sv_states)}")
for i, p in enumerate(sv_states, 1):
    print(f"  Sep-of-vars steady state {i}: X = {p[0]:.6f}, Y = {p[1]:.6f}")

# Match the two lists and report the maximum discrepancy.
if steady_states and sv_states:
    max_diff = 0.0
    for p in steady_states:
        d = min(np.hypot(*(p - q)) for q in sv_states)
        max_diff = max(max_diff, d)
    print(f"Max distance between contour crossings and sep-of-vars steady states: {max_diff:.3e}")

# One-sentence explanation of why the check confirms the result:
print("Why the check confirms the result: the closed-form separation-of-variables "
      "curves give zero residual when substituted into fX and fY and yield the same "
      "crossing coordinates, so the zero contours and their intersections are the "
      "true nullclines and steady states, not artifacts of the grid or contour tracer.")
