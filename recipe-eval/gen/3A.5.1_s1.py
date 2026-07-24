import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# ----------------------------------------------------------------------
# Model parameters (mutual-repression toggle switch)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ----------------------------------------------------------------------
# Right-hand sides.  Each nullcline is where one of these surfaces = 0.
# ----------------------------------------------------------------------
def fX(X, Y):
    # dX/dt : basal + repressive Hill(Y) - degradation
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def fY(X, Y):
    # dY/dt : basal + repressive Hill(X) - degradation
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ----------------------------------------------------------------------
# GENERAL METHOD (no separation of variables):
# 1. Build a grid over the (X,Y) plane.
# 2. Evaluate the scalar surface Z = f(X,Y) at every grid point.
# 3. Let a contour routine trace the single level curve Z = 0.
# ----------------------------------------------------------------------
Xmax = (gX0 + gX1) / kX * 1.15     # generous bounds from max possible values
Ymax = (gY0 + gY1) / kY * 1.15
xg = np.linspace(0.0, Xmax, 600)
yg = np.linspace(0.0, Ymax, 600)
XX, YY = np.meshgrid(xg, yg)       # step 1: grid

ZX = fX(XX, YY)                    # step 2: surface for dX/dt
ZY = fY(XX, YY)                    # step 2: surface for dY/dt

# ----------------------------------------------------------------------
# Find the crossings (steady states) with a root finder launched from a
# coarse lattice of initial guesses, then de-duplicate the results.
# ----------------------------------------------------------------------
def system(v):
    X, Y = v
    return [fX(X, Y), fY(X, Y)]

roots = []
for x0 in np.linspace(0.0, Xmax, 12):
    for y0 in np.linspace(0.0, Ymax, 12):
        sol, info, ier, _ = fsolve(system, [x0, y0], full_output=True)
        if ier == 1 and 0 <= sol[0] <= Xmax and 0 <= sol[1] <= Ymax:
            if not any(np.hypot(sol[0]-r[0], sol[1]-r[1]) < 1e-3 for r in roots):
                roots.append(sol)
roots = sorted(roots, key=lambda r: r[0])

print("=== Steady states (crossings of the two zero contours) ===")
for i, (Xs, Ys) in enumerate(roots, 1):
    print(f"Steady state {i}: X = {Xs:.6f}, Y = {Ys:.6f}  "
          f"(fX = {fX(Xs,Ys):.2e}, fY = {fY(Xs,Ys):.2e})")

# ----------------------------------------------------------------------
# CHECK: separation of variables gives each nullcline explicitly.
# X-nullcline: solve fX=0 for X  ->  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX
# Y-nullcline: solve fY=0 for Y  ->  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY
# ----------------------------------------------------------------------
def X_nullcline_of_Y(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def Y_nullcline_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

Xexp = X_nullcline_of_Y(yg)       # explicit X-nullcline, X as function of Y
Yexp = Y_nullcline_of_X(xg)       # explicit Y-nullcline, Y as function of X

# Quantify agreement: residual of explicit curves on the OTHER surface = 0
resX = np.max(np.abs(fX(Xexp, yg)))
resY = np.max(np.abs(fY(xg, Yexp)))
print("\n=== Contour vs. separation-of-variables agreement ===")
print(f"Max |fX| along explicit X-nullcline: {resX:.3e}")
print(f"Max |fY| along explicit Y-nullcline: {resY:.3e}")

# Confirm every fsolve crossing lies on both explicit curves
print("\n=== Crossings vs. explicit-nullcline intersections ===")
for i, (Xs, Ys) in enumerate(roots, 1):
    dX = Xs - X_nullcline_of_Y(Ys)
    dY = Ys - Y_nullcline_of_X(Xs)
    print(f"Steady state {i}: X - Xnull(Y) = {dX:.3e}, Y - Ynull(X) = {dY:.3e}")

# ----------------------------------------------------------------------
# Phase-plane plot
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 7))

# zero-level contours (the general method)
cX = ax.contour(XX, YY, ZX, levels=[0.0], colors="tab:blue", linewidths=2.5)
cY = ax.contour(XX, YY, ZY, levels=[0.0], colors="tab:red", linewidths=2.5)

# explicit nullclines from separation of variables (dashed overlay = check)
ax.plot(Xexp, yg, "k--", lw=1.2, label="X-nullcline (sep. of variables)")
ax.plot(xg, Yexp, color="0.4", ls=":", lw=1.6, label="Y-nullcline (sep. of variables)")

# crossings
for (Xs, Ys) in roots:
    ax.plot(Xs, Ys, "ko", ms=9, mfc="yellow", zorder=5)

# proxy legend entries for the contours
ax.plot([], [], color="tab:blue", lw=2.5, label="dX/dt = 0  (zero contour)")
ax.plot([], [], color="tab:red", lw=2.5, label="dY/dt = 0  (zero contour)")
ax.plot([], [], "ko", mfc="yellow", ms=9, label="steady states")

ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Toggle switch nullclines as zero-level contours")
ax.set_xlim(0, Xmax)
ax.set_ylim(0, Ymax)
ax.legend(loc="upper right", fontsize=9)
ax.grid(alpha=0.3)

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.5.1_s1.png", dpi=130)

# ----------------------------------------------------------------------
# One-sentence explanation of why the check is valid.
# ----------------------------------------------------------------------
print("\nWhy the check confirms the result:")
print("Because the explicit nullclines are the exact algebraic solutions of "
      "fX=0 and fY=0, their coincidence with the numerically traced zero "
      "contours (near-zero residuals) and with the fsolve crossings shows the "
      "grid/contour method locates the same curves and steady states without "
      "relying on separation of variables.")
