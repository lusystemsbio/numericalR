import numpy as np
from scipy.optimize import fsolve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Model parameters (toggle switch: X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10   # dX/dt parameters
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12   # dY/dt parameters

# Right-hand sides.  fX = dX/dt, fY = dY/dt.
# Basal rate + repressive Hill function of the OTHER gene - linear degradation.
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def fY(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ----------------------------------------------------------------------
# 1) GENERAL METHOD: evaluate each RHS as a surface Z on a grid,
#    then trace its zero-level contour.  No separation of variables.
# ----------------------------------------------------------------------
xmax, ymax = 700.0, 500.0
xs = np.linspace(0.0, xmax, 400)
ys = np.linspace(0.0, ymax, 400)
XX, YY = np.meshgrid(xs, ys)          # grid over the (X, Y) plane
ZX = fX(XX, YY)                       # surface for the X-nullcline
ZY = fY(XX, YY)                       # surface for the Y-nullcline

fig, ax = plt.subplots(figsize=(8, 6))
# The nullcline is exactly the set where the surface passes through zero,
# so we let the contour routine trace the level-0 isoline of each surface.
cX = ax.contour(XX, YY, ZX, levels=[0.0], colors="crimson", linewidths=2)
cY = ax.contour(XX, YY, ZY, levels=[0.0], colors="royalblue", linewidths=2)

# ----------------------------------------------------------------------
# 2) Crossings of the two zero contours = steady states of the system.
#    Found generically by solving (fX, fY) = (0, 0) from a grid of guesses.
# ----------------------------------------------------------------------
def rhs(p):
    X, Y = p
    return [fX(X, Y), fY(X, Y)]

roots = []
for gx in np.linspace(20, xmax, 8):
    for gy in np.linspace(20, ymax, 8):
        sol, info, ier, msg = fsolve(rhs, [gx, gy], full_output=True)
        if ier == 1 and 0 <= sol[0] <= xmax and 0 <= sol[1] <= ymax:
            # keep only genuinely new roots (dedupe)
            if not any(np.hypot(sol[0]-r[0], sol[1]-r[1]) < 1e-3 for r in roots):
                roots.append((sol[0], sol[1]))
roots.sort()

ax.plot([r[0] for r in roots], [r[1] for r in roots],
        "ko", ms=9, mfc="yellow", mec="k", label="crossings (steady states)")

# ----------------------------------------------------------------------
# 3) CHECK: separation of variables gives each nullcline explicitly.
#    dX/dt = 0  ->  X = [gX0 + gX1/(1+(Y/Yth)^nY)] / kX   (X as function of Y)
#    dY/dt = 0  ->  Y = [gY0 + gY1/(1+(X/Xth)^nX)] / kY   (Y as function of X)
# ----------------------------------------------------------------------
def X_null_of_Y(Y):   # solved X-nullcline
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX
def Y_null_of_X(X):   # solved Y-nullcline
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

ax.plot(X_null_of_Y(ys), ys, "k--", lw=1, label="X-nullcline (sep. of vars)")
ax.plot(xs, Y_null_of_X(xs), "k:", lw=1, label="Y-nullcline (sep. of vars)")

ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Toggle switch nullclines as zero contours, with crossings")
ax.set_xlim(0, xmax); ax.set_ylim(0, ymax)
# proxy legend entries for the contour colors
ax.plot([], [], color="crimson", lw=2, label="X-nullcline (zero contour)")
ax.plot([], [], color="royalblue", lw=2, label="Y-nullcline (zero contour)")
ax.legend(loc="upper right", fontsize=8)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.5.1_s2.png")

# ----------------------------------------------------------------------
# 4) Quantitative confirmation that the two methods agree.
# ----------------------------------------------------------------------
# (a) The zero-contour vertices should lie on the separation-of-variables curves.
def contour_max_residual(cset, on_curve, use_X):
    # residual = signed distance of contour vertices from the explicit curve
    res = 0.0
    for path in cset.allsegs[0]:
        if len(path) == 0:
            continue
        Xv, Yv = path[:, 0], path[:, 1]
        pred = on_curve(Yv) if use_X else on_curve(Xv)
        target = Xv if use_X else Yv
        res = max(res, np.max(np.abs(target - pred)))
    return res

resX = contour_max_residual(cX, X_null_of_Y, use_X=True)
resY = contour_max_residual(cY, Y_null_of_X, use_X=False)

print(f"Max |X_contour - X_null(Y)| over X-nullcline vertices: {resX:.4e}")
print(f"Max |Y_contour - Y_null(X)| over Y-nullcline vertices: {resY:.4e}")

# (b) Steady states via substitution (separation of variables): find X where
#     X = X_null_of_Y( Y_null_of_X(X) ), i.e. self-consistency of both curves.
def sub_eq(X):
    return X - X_null_of_Y(Y_null_of_X(X))

sub_roots = []
for gx in np.linspace(20, xmax, 200):
    Xr = fsolve(sub_eq, gx)[0]
    if 0 <= Xr <= xmax and not any(abs(Xr - r) < 1e-3 for r in sub_roots):
        sub_roots.append(Xr)
sub_roots.sort()
sub_states = [(Xr, Y_null_of_X(Xr)) for Xr in sub_roots]

print("\nSteady states from zero-contour crossings (fsolve on fX=fY=0):")
for X, Y in roots:
    print(f"  X = {X:10.4f}   Y = {Y:10.4f}   fX = {fX(X,Y):.2e}   fY = {fY(X,Y):.2e}")

print("\nSteady states from separation of variables (substitution):")
for X, Y in sub_states:
    print(f"  X = {X:10.4f}   Y = {Y:10.4f}")

print("\nMatch between the two sets of steady states (max coordinate difference):")
for (Xa, Ya) in roots:
    d = min(np.hypot(Xa - Xb, Ya - Yb) for (Xb, Yb) in sub_states)
    print(f"  crossing (X={Xa:.4f}, Y={Ya:.4f}) -> nearest sep-of-vars state distance = {d:.4e}")

# One-sentence explanation of why the check confirms the result:
print("\nWhy the check confirms the result:")
print("Because the zero-level contour of fX (resp. fY) is by definition the set "
      "where dX/dt=0 (resp. dY/dt=0), it must coincide point-for-point with the "
      "curve obtained by algebraically solving that same equation, and their "
      "intersections are exactly the points where both derivatives vanish, i.e. "
      "the steady states.")
