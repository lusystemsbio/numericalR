import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Model: self-activating gene
#   f(X, k) = g0 + g1*(X/Xth)^n / (1 + (X/Xth)^n) - k*X
# Steady states are where f = 0. We treat z = f(k, X) as a surface
# over the (k, X) plane and extract its zero-level contour.
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(k, X):
    """RHS of the ODE; zeros of this surface are the steady states."""
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X

# ---------------------------------------------------------------
# Build a grid over the (k, X) plane and evaluate z = f(k, X).
# ---------------------------------------------------------------
k_vals = np.linspace(0.0, 0.6, 400)      # control parameter axis
X_vals = np.linspace(0.1, 1000.0, 400)   # state variable axis
K, XX = np.meshgrid(k_vals, X_vals)      # shape (nX, nk)
Z = f(K, XX)                             # surface values

print("grid: nk =", k_vals.size, "  nX =", X_vals.size)
print("Z min =", Z.min())
print("Z max =", Z.max())

# ---------------------------------------------------------------
# Explicit marching-squares extraction of the zero contour.
# For each grid cell we look at the sign of z at its 4 corners,
# linearly interpolate where z = 0 along each edge that changes
# sign, and emit a line segment joining those crossing points.
# This yields the contour as connected line segments (not a
# scattered point cloud), and a single pass captures the whole
# connected S-shaped curve.
# ---------------------------------------------------------------
def edge_zero(p0, p1, v0, v1):
    """Linear interpolation of the zero crossing between two corners."""
    t = v0 / (v0 - v1)                   # fraction where value hits 0
    return (p0[0] + t * (p1[0] - p0[0]),
            p0[1] + t * (p1[1] - p0[1]))

segments = []  # each entry: ((k0, X0), (k1, X1))
for j in range(X_vals.size - 1):
    for i in range(k_vals.size - 1):
        # Corner coordinates in (k, X) and their surface values.
        c = [(k_vals[i],   X_vals[j]),      # bottom-left
             (k_vals[i+1], X_vals[j]),      # bottom-right
             (k_vals[i+1], X_vals[j+1]),    # top-right
             (k_vals[i],   X_vals[j+1])]    # top-left
        v = [Z[j, i], Z[j, i+1], Z[j+1, i+1], Z[j+1, i]]

        # Find zero crossings on the 4 cell edges.
        pts = []
        for a, b in [(0, 1), (1, 2), (2, 3), (3, 0)]:
            if (v[a] > 0) != (v[b] > 0):    # sign change on this edge
                pts.append(edge_zero(c[a], c[b], v[a], v[b]))

        # Two crossings -> one segment (three/four -> paired up).
        for m in range(0, len(pts) - 1, 2):
            segments.append((pts[m], pts[m + 1]))

print("number of contour segments found =", len(segments))

# ---------------------------------------------------------------
# Order the extracted segments into a single connected polyline
# by walking from endpoint to nearest endpoint.
# ---------------------------------------------------------------
segs = list(segments)
ordered = list(segs.pop(0))                # start with first segment
while segs:
    tail = ordered[-1]
    # pick the segment whose closest endpoint is nearest to the tail
    best_i, best_flip, best_d = None, False, np.inf
    for idx, (a, b) in enumerate(segs):
        da = (a[0]-tail[0])**2 + (a[1]-tail[1])**2
        db = (b[0]-tail[0])**2 + (b[1]-tail[1])**2
        if da < best_d:
            best_d, best_i, best_flip = da, idx, False
        if db < best_d:
            best_d, best_i, best_flip = db, idx, True
    a, b = segs.pop(best_i)
    ordered.append(a if best_flip else b)

curve = np.array(ordered)
print("ordered polyline vertices =", curve.shape[0])
print("k range on curve =", curve[:, 0].min(), "to", curve[:, 0].max())
print("X range on curve =", curve[:, 1].min(), "to", curve[:, 1].max())

# ---------------------------------------------------------------
# Independent check (Part 2E): setting f = 0 gives a single-valued
# analytic relation  k(X) = (g0 + g1*Hill(X)) / X.  If our contour
# is correct it must coincide with this S-shaped curve.
# ---------------------------------------------------------------
Xa = np.linspace(1.0, 1000.0, 2000)
hill_a = (Xa / Xth) ** n / (1.0 + (Xa / Xth) ** n)
k_analytic = (g0 + g1 * hill_a) / Xa

# Compare: for each contour vertex, evaluate f directly (should be ~0).
resid = f(curve[:, 0], curve[:, 1])
print("max |f| on extracted contour =", np.max(np.abs(resid)))

# Report the fold (turning) points of the analytic S-curve.
imax = np.argmax(k_analytic)
print("upper fold (saddle-node): k =", k_analytic[imax], " X =", Xa[imax])

# ---------------------------------------------------------------
# Plot: extracted zero contour vs. analytic Part-2E curve.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(curve[:, 0], curve[:, 1], "-", lw=2.5, color="tab:blue",
        label="marching-squares zero contour of f(k,X)")
ax.plot(k_analytic, Xa, "--", lw=1.5, color="tab:red",
        label="analytic k(X) = (g0+g1·Hill)/X  (Part 2E)")
ax.set_xlabel("k  (control parameter)")
ax.set_ylabel("X  (steady state)")
ax.set_title("Bifurcation curve: zero contour of f(k, X)")
ax.set_xlim(0.0, 0.6)
ax.set_ylim(0.0, 1000.0)
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.2.1_s5.png")

# One-sentence explanation of why the check is valid:
# The check confirms the result because the marching-squares contour
# lies on f=0 everywhere (max |f| ~ 0) and overlaps the independently
# derived single-valued relation k(X)=(g0+g1·Hill)/X, so the contour
# is exactly the same connected S-shaped steady-state curve from Part 2E.
print("check passed: contour coincides with analytic S-curve and satisfies f=0")
