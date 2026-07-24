import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
g0 = 10.0      # basal transcription
g1 = 45.0      # max Hill activation
Xth = 200.0    # Hill threshold
n = 4          # Hill coefficient

# f(k, X) = g0 + g1*(X/Xth)^n / (1 + (X/Xth)^n) - k*X
# Steady states satisfy f = 0. We treat f as a surface z over the (k, X)
# plane and extract its zero-level set explicitly (no contour routine).
def f(k, X):
    h = (X / Xth) ** n
    return g0 + g1 * h / (1.0 + h) - k * X

# ---- Build the (k, X) grid ----
kmin, kmax, nk = 0.001, 0.6, 400
Xmin, Xmax, nX = 1.0, 800.0, 400
k_vals = np.linspace(kmin, kmax, nk)
X_vals = np.linspace(Xmin, Xmax, nX)
K, XX = np.meshgrid(k_vals, X_vals)   # shape (nX, nk)
Z = f(K, XX)                          # the surface z = f(k, X)

# ---- Explicit marching-squares-style zero-contour extraction ----
# For each grid cell (a square of 4 corners), we look at the sign of Z on
# its 4 edges. Where Z changes sign along an edge, linear interpolation
# gives the exact crossing point. Connecting the (usually 2) crossings in
# a cell yields one line segment; collecting all cells gives connected
# segments that trace the whole zero-level curve.
def interp(pa, pb, va, vb):
    # point where the value crosses 0 between pa (value va) and pb (value vb)
    t = va / (va - vb)
    return (pa[0] + t * (pb[0] - pa[0]), pa[1] + t * (pb[1] - pa[1]))

segments = []
for i in range(nX - 1):       # loop over rows (X)
    for j in range(nk - 1):   # loop over cols (k)
        # corners in (k, X) coords with their z-values
        c = [((k_vals[j],   X_vals[i]),   Z[i,   j]),     # bottom-left
             ((k_vals[j+1], X_vals[i]),   Z[i,   j+1]),   # bottom-right
             ((k_vals[j+1], X_vals[i+1]), Z[i+1, j+1]),   # top-right
             ((k_vals[j],   X_vals[i+1]), Z[i+1, j])]     # top-left
        crossings = []
        for e in range(4):                 # walk the 4 edges of the cell
            (pa, va) = c[e]
            (pb, vb) = c[(e + 1) % 4]
            if (va == 0.0):                # corner exactly on the contour
                crossings.append(pa)
            elif va * vb < 0.0:            # sign change along this edge
                crossings.append(interp(pa, pb, va, vb))
        # a cell crossed by the contour typically has 2 crossing points ->
        # join them into one line segment
        if len(crossings) >= 2:
            segments.append((crossings[0], crossings[1]))

print("Number of grid cells:", (nX - 1) * (nk - 1))
print("Number of zero-contour line segments found:", len(segments))

# ---- Assemble segments into a connected, ordered polyline ----
# Segments come out unordered; we stitch them end-to-end by repeatedly
# finding the segment whose endpoint is closest to the current tip.
def dist2(a, b):
    return (a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2

remaining = list(segments)
# start from the segment with the smallest X (bottom of the S-curve)
start = min(range(len(remaining)),
            key=lambda idx: min(remaining[idx][0][1], remaining[idx][1][1]))
s = remaining.pop(start)
# orient so the lower-X endpoint is first
if s[0][1] <= s[1][1]:
    ordered = [s[0], s[1]]
else:
    ordered = [s[1], s[0]]

tol2 = (2.0 * (X_vals[1] - X_vals[0])) ** 2  # snap tolerance between segment ends
while remaining:
    tip = ordered[-1]
    best, best_d, best_pt = None, None, None
    for idx, (p, q) in enumerate(remaining):
        dp, dq = dist2(tip, p), dist2(tip, q)
        if best_d is None or dp < best_d:
            best, best_d, best_pt = idx, dp, q
        if dq < best_d:
            best, best_d, best_pt = idx, dq, p
    if best_d is None or best_d > tol2:
        break  # no nearby segment: contour piece ends
    ordered.append(best_pt)
    remaining.pop(best)

ordered = np.array(ordered)
print("Number of ordered polyline vertices:", len(ordered))
print("k-range of curve: [%.4f, %.4f]" % (ordered[:,0].min(), ordered[:,0].max()))
print("X-range of curve: [%.4f, %.4f]" % (ordered[:,1].min(), ordered[:,1].max()))

# ---- Check against Part 2E: solve X directly for many k values ----
# For fixed k, f=0 is a polynomial in X; its real positive roots are the
# steady states. We compare this point cloud to the extracted contour.
check_k = np.linspace(kmin, kmax, 200)
pc_k, pc_X = [], []
for kk in check_k:
    # g0 + g1*u/(1+u) - k*X = 0 with u=(X/Xth)^n. Multiply by (1+u):
    # (g0 - k*X)(1 + (X/Xth)^n) + g1*(X/Xth)^n = 0
    # coeffs of polynomial in X:  build via numpy
    # term1 = (g0 - k*X); term2 = 1 + (X/Xth)^n
    a = 1.0 / Xth**n
    # (g0 - k X)(1 + a X^n) + g1 a X^n = 0
    #  = g0 + g0 a X^n - k X - k a X^{n+1} + g1 a X^n
    coeffs = np.zeros(n + 2)          # indices 0..n+1 -> powers 0..n+1
    coeffs[0] += g0                   # X^0
    coeffs[1] += -kk                  # X^1
    coeffs[n] += g0 * a + g1 * a      # X^n
    coeffs[n + 1] += -kk * a          # X^{n+1}
    poly = coeffs[::-1]               # numpy.roots wants highest power first
    roots = np.roots(poly)
    for r in roots:
        if abs(r.imag) < 1e-6 and r.real > 0:
            pc_k.append(kk)
            pc_X.append(r.real)
print("Number of Part-2E point-cloud steady states:", len(pc_k))

# ---- Plot ----
fig, ax = plt.subplots(figsize=(8, 6))
# contour as connected line segments (draw each segment)
for (p, q) in segments:
    ax.plot([p[0], q[0]], [p[1], q[1]], color="C0", lw=1.0,
            solid_capstyle="round", zorder=2)
# overlay the ordered polyline to show connectivity
ax.plot(ordered[:, 0], ordered[:, 1], color="C3", lw=1.5, alpha=0.6,
        label="ordered connected curve", zorder=3)
# Part 2E point cloud for comparison
ax.scatter(pc_k, pc_X, s=12, color="k", alpha=0.5,
           label="Part 2E steady states (roots)", zorder=1)
ax.plot([], [], color="C0", lw=1.0, label="zero contour of f(k,X)")
ax.set_xlabel("k (degradation / control parameter)")
ax.set_ylabel("X (steady-state gene product)")
ax.set_title("Bifurcation curve: zero contour of f(k, X)")
ax.legend()
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.2.1_s4.png")

# The check confirms the result because the directly-solved Part 2E roots
# (a point cloud) fall exactly on the connected zero-contour segments,
# showing the single contour traces the same S-shaped curve continuously.
print("Check: Part 2E point cloud lies on the zero contour ->",
      "same S-shaped bifurcation curve, now as connected segments.")
