import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Right-hand sides.  A nullcline is the set where one of these is zero.
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X   # dX/dt

def fY(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y   # dY/dt

# ----------------------------------------------------------------------
# GENERAL METHOD (no separation of variables):
# evaluate the surface Z = f(X,Y) on a grid and trace the Z = 0 contour.
# ----------------------------------------------------------------------
Xg = np.linspace(0.0, 600.0, 500)      # X axis samples
Yg = np.linspace(0.0, 450.0, 500)      # Y axis samples
XX, YY = np.meshgrid(Xg, Yg)           # XX[j,i]=Xg[i], YY[j,i]=Yg[j]
ZX = fX(XX, YY)                        # surface for the X-nullcline
ZY = fY(XX, YY)                        # surface for the Y-nullcline

# Explicit marching-squares: for every grid cell, find where the surface
# changes sign along the four edges, linearly interpolate the exact zero
# location, and join the two crossings into a short segment of the contour.
def zero_contour_segments(Z, X, Y):
    segs = []
    ny, nx = Z.shape
    for j in range(ny - 1):
        for i in range(nx - 1):
            # corner values and coordinates of this cell
            v00, v10 = Z[j, i],     Z[j, i + 1]      # bottom-left, bottom-right
            v01, v11 = Z[j + 1, i], Z[j + 1, i + 1]  # top-left,   top-right
            x0, x1 = X[i], X[i + 1]
            y0, y1 = Y[j], Y[j + 1]
            pts = []
            # bottom edge: sign change between v00 and v10
            if (v00 > 0) != (v10 > 0):
                t = v00 / (v00 - v10)
                pts.append((x0 + t * (x1 - x0), y0))
            # top edge
            if (v01 > 0) != (v11 > 0):
                t = v01 / (v01 - v11)
                pts.append((x0 + t * (x1 - x0), y1))
            # left edge
            if (v00 > 0) != (v01 > 0):
                t = v00 / (v00 - v01)
                pts.append((x0, y0 + t * (y1 - y0)))
            # right edge
            if (v10 > 0) != (v11 > 0):
                t = v10 / (v10 - v11)
                pts.append((x1, y0 + t * (y1 - y0)))
            # 2 crossings -> one segment; 4 -> ambiguous saddle, split in pairs
            if len(pts) == 2:
                segs.append((pts[0], pts[1]))
            elif len(pts) == 4:
                segs.append((pts[0], pts[1]))
                segs.append((pts[2], pts[3]))
    return segs

segX = zero_contour_segments(ZX, Xg, Yg)   # traced X-nullcline (fX=0)
segY = zero_contour_segments(ZY, Xg, Yg)   # traced Y-nullcline (fY=0)

# ----------------------------------------------------------------------
# STEADY STATES = nullcline crossings.
# Using separation of variables the nullclines are explicit functions:
#   fX=0  ->  X = FX(Y) = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX
#   fY=0  ->  Y = FY(X) = (gY0 + gY1/(1+(X/Xth)^nX)) / kY
# A fixed point satisfies X = FX(FY(X)); find its roots by scanning + bisection.
# ----------------------------------------------------------------------
def FX(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def FY(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

def g(X):
    return FX(FY(X)) - X

Xscan = np.linspace(1.0, 600.0, 4000)
gvals = g(Xscan)
roots = []
for i in range(len(Xscan) - 1):
    if gvals[i] == 0.0:
        roots.append(Xscan[i])
    elif (gvals[i] > 0) != (gvals[i + 1] > 0):    # sign change -> bracketed root
        a, b = Xscan[i], Xscan[i + 1]
        for _ in range(80):                        # bisection refinement
            m = 0.5 * (a + b)
            if (g(a) > 0) != (g(m) > 0):
                b = m
            else:
                a = m
        roots.append(0.5 * (a + b))

steady = [(Xs, FY(Xs)) for Xs in roots]

# ----------------------------------------------------------------------
# Numerical checks
# ----------------------------------------------------------------------
print("=== Steady states (nullcline crossings) ===")
for k, (Xs, Ys) in enumerate(steady):
    print(f"steady_state_{k+1}_X = {Xs:.6f}")
    print(f"steady_state_{k+1}_Y = {Ys:.6f}")
    print(f"steady_state_{k+1}_fX_residual = {fX(Xs, Ys):.3e}")
    print(f"steady_state_{k+1}_fY_residual = {fY(Xs, Ys):.3e}")

# Residual of fX along the traced X-nullcline segment midpoints (should be ~0),
# and residual of the separation-of-variables curve, to confirm they agree.
midsX = np.array([[(p[0] + q[0]) / 2, (p[1] + q[1]) / 2] for p, q in segX])
midsY = np.array([[(p[0] + q[0]) / 2, (p[1] + q[1]) / 2] for p, q in segY])
print("max_|fX|_on_traced_X_nullcline =",
      f"{np.max(np.abs(fX(midsX[:,0], midsX[:,1]))):.3e}")
print("max_|fY|_on_traced_Y_nullcline =",
      f"{np.max(np.abs(fY(midsY[:,0], midsY[:,1]))):.3e}")

# Compare traced contour to the separation-of-variables curve for the X-nullcline:
# for each traced midpoint, X should equal FX(Y).
devX = np.max(np.abs(midsX[:, 0] - FX(midsX[:, 1])))
devY = np.max(np.abs(midsY[:, 1] - FY(midsY[:, 0])))
print("max_dev_X_contour_vs_FX(Y) =", f"{devX:.3e}")
print("max_dev_Y_contour_vs_FY(X) =", f"{devY:.3e}")

# ----------------------------------------------------------------------
# Plots: (left) zero-contour method, (right) check vs separation of variables
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

# Left: nullclines drawn purely as traced zero contours, plus their crossings.
ax1.add_collection(LineCollection(segX, colors="tab:blue", linewidths=1.5,
                                  label="X-nullcline (fX=0)"))
ax1.add_collection(LineCollection(segY, colors="tab:red", linewidths=1.5,
                                  label="Y-nullcline (fY=0)"))
ax1.plot([], [], color="tab:blue", label="X-nullcline (fX=0)")   # legend proxies
ax1.plot([], [], color="tab:red", label="Y-nullcline (fY=0)")
for Xs, Ys in steady:
    ax1.plot(Xs, Ys, "ko", ms=8, zorder=5)
ax1.set_title("Nullclines as zero-level contours (general method)")
ax1.set_xlabel("X"); ax1.set_ylabel("Y")
ax1.set_xlim(0, 600); ax1.set_ylim(0, 450); ax1.legend(); ax1.grid(alpha=0.3)

# Right: overlay separation-of-variables curves on the traced contours.
Yc = np.linspace(0, 450, 400)
Xc = np.linspace(0, 600, 400)
ax2.add_collection(LineCollection(segX, colors="tab:blue", linewidths=3.0, alpha=0.35))
ax2.add_collection(LineCollection(segY, colors="tab:red", linewidths=3.0, alpha=0.35))
ax2.plot(FX(Yc), Yc, "b--", lw=1.2, label="X-nullcline: X=FX(Y) (sep. of vars)")
ax2.plot(Xc, FY(Xc), "r--", lw=1.2, label="Y-nullcline: Y=FY(X) (sep. of vars)")
for Xs, Ys in steady:
    ax2.plot(Xs, Ys, "ko", ms=8, zorder=5)
ax2.set_title("Check: zero contours (thick) vs separation of variables (dashed)")
ax2.set_xlabel("X"); ax2.set_ylabel("Y")
ax2.set_xlim(0, 600); ax2.set_ylim(0, 450); ax2.legend(); ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.5.1_s3.png")

# One-sentence explanation of why the check is valid:
print("EXPLANATION: The zero contours and the separation-of-variables curves are "
      "two independent renderings of the exact same equations fX=0 and fY=0, so "
      "their coincidence (tiny deviations) and identical crossing points confirm "
      "that the contour method traced the true nullclines and steady states.")
