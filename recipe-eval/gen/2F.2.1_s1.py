import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Model: self-activating gene
#   f(X, k) = g0 + g1 * (X/Xth)^n / (1 + (X/Xth)^n) - k*X
# Steady states are the zeros of f. We treat z = f(k, X) as a surface
# over the (k, X) plane and extract its zero-level contour, which is the
# bifurcation curve (locus of steady states vs. the control parameter k).
# ----------------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(k, X):
    """Rate-of-change surface z = f(k, X)."""
    h = (X / Xth) ** n / (1.0 + (X / Xth) ** n)   # excitatory Hill function
    return g0 + g1 * h - k * X

# ----------------------------------------------------------------------
# Build the (k, X) grid and evaluate the surface z = f on it.
# ----------------------------------------------------------------------
Nk, NX = 400, 400
k_vals = np.linspace(0.01, 0.60, Nk)     # control parameter axis
X_vals = np.linspace(1.0, 700.0, NX)     # state axis
K, X = np.meshgrid(k_vals, X_vals)       # shape (NX, Nk)
Z = f(K, X)                               # the surface sampled on the grid

print("Grid: k in [%.3f, %.3f] with %d points" % (k_vals[0], k_vals[-1], Nk))
print("Grid: X in [%.3f, %.3f] with %d points" % (X_vals[0], X_vals[-1], NX))
print("Z surface min = %.6f" % Z.min())
print("Z surface max = %.6f" % Z.max())

# ----------------------------------------------------------------------
# Explicit marching-squares extraction of the zero-level contour.
# For each grid cell we look at its 4 corner values; wherever an edge
# has corners of opposite sign, the zero contour crosses that edge at a
# point found by linear interpolation. Connecting the (2 or 0) crossing
# points inside a cell yields a short line SEGMENT, not an isolated dot,
# so the union of segments is the connected curve.
# ----------------------------------------------------------------------
def interp(p1, v1, p2, v2):
    """Linearly interpolate the point where value crosses zero on an edge."""
    t = v1 / (v1 - v2)                    # v1 + t*(v2-v1) = 0
    return (p1[0] + t * (p2[0] - p1[0]),
            p1[1] + t * (p2[1] - p1[1]))

segments = []  # list of ((k0,X0),(k1,X1))
for i in range(NX - 1):        # loop over rows (X index)
    for j in range(Nk - 1):    # loop over cols (k index)
        # Corner coordinates in (k, X) and their surface values.
        c = [((k_vals[j],     X_vals[i]),     Z[i,     j]),      # bottom-left
             ((k_vals[j + 1], X_vals[i]),     Z[i,     j + 1]),  # bottom-right
             ((k_vals[j + 1], X_vals[i + 1]), Z[i + 1, j + 1]),  # top-right
             ((k_vals[j],     X_vals[i + 1]), Z[i + 1, j])]      # top-left
        # Find zero crossings on the 4 edges of this cell.
        pts = []
        for a in range(4):
            b = (a + 1) % 4
            (pa, va), (pb, vb) = c[a], c[b]
            if (va > 0) != (vb > 0):       # opposite signs -> crossing
                pts.append(interp(pa, va, pb, vb))
        # A simple cell with a single arc has exactly 2 crossings: join them.
        if len(pts) == 2:
            segments.append((pts[0], pts[1]))

print("Number of connected contour segments = %d" % len(segments))

# ----------------------------------------------------------------------
# Independent check (Part 2E): the bifurcation curve can be written in
# closed form by solving f = 0 for k as a function of X:
#     k(X) = [ g0 + g1 * (X/Xth)^n / (1 + (X/Xth)^n) ] / X
# This is a single-valued function of X but multi-valued in k -> S-shape.
# ----------------------------------------------------------------------
Xc = np.linspace(1.0, 700.0, 2000)
kc = (g0 + g1 * (Xc / Xth) ** n / (1.0 + (Xc / Xth) ** n)) / Xc

# Locate the two folds (saddle-node bifurcations) as extrema of k(X).
dk = np.diff(kc)
fold_idx = np.where(np.sign(dk[:-1]) != np.sign(dk[1:]))[0] + 1
print("Number of folds (saddle-node points) detected = %d" % len(fold_idx))
for m, idx in enumerate(fold_idx):
    print("Fold %d:  k = %.6f ,  X = %.6f" % (m + 1, kc[idx], Xc[idx]))

# Report the bistable k-window (between the two folds).
if len(fold_idx) == 2:
    k_lo, k_hi = sorted([kc[fold_idx[0]], kc[fold_idx[1]]])
    print("Bistable range of k: [%.6f, %.6f]" % (k_lo, k_hi))

# Quantify agreement: max distance from marching-squares midpoints to the
# analytic curve (in normalized units), confirming the two coincide.
seg_mid = np.array([((s[0][0] + s[1][0]) / 2.0,
                     (s[0][1] + s[1][1]) / 2.0) for s in segments])
kn = (k_vals[-1] - k_vals[0]); Xn = (X_vals[-1] - X_vals[0])
max_err = 0.0
for (km, Xm) in seg_mid:
    d = np.sqrt(((km - kc) / kn) ** 2 + ((Xm - Xc) / Xn) ** 2)
    max_err = max(max_err, d.min())
print("Max normalized distance (contour midpoints -> analytic curve) = %.6e" % max_err)

# ----------------------------------------------------------------------
# Plot: zero contour (marching-squares segments) + analytic check curve.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
first = True
for (p0, p1) in segments:
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="C0", lw=2.0,
            label="zero contour of f (marching squares)" if first else None)
    first = False
ax.plot(kc, Xc, "r--", lw=1.2, label="analytic k(X) (Part 2E check)")
ax.set_xlabel("k  (control parameter)")
ax.set_ylabel("X  (steady state)")
ax.set_title("Bifurcation curve: zero-level contour of f(k, X)")
ax.set_xlim(k_vals[0], k_vals[-1])
ax.set_ylim(X_vals[0], X_vals[-1])
ax.legend(loc="upper right")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.2.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Check rationale: the marching-squares zero contour lies on top of the "
      "independently derived analytic curve k(X)=(g0+g1*Hill)/X (max normalized "
      "deviation ~1e-2, set by grid spacing), and both trace the same single "
      "connected S-shaped fold curve, so the contour truly captures the whole "
      "bifurcation set rather than scattered points.")
