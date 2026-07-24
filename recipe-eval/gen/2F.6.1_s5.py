import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Model: self-activating gene
#   f(X,k) = g0 + g1*Hill(X) - k*X ,  Hill(X) = u/(1+u), u=(X/Xth)^n
# Equilibria are f(X,k)=0. We trace this curve in the (k, X) plane.
# ----------------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def hill(X):
    u = (X / Xth) ** n
    return u / (1.0 + u)

def f(X, k):
    return g0 + g1 * hill(X) - k * X

def fX(X, k):
    # partial derivative wrt X.  dHill/dX = (n*u/X)/(1+u)^2
    u = (X / Xth) ** n
    dHill = (n * u / X) / (1.0 + u) ** 2
    return g1 * dHill - k

def fk(X, k):
    return -X                      # partial derivative wrt k

# ----------------------------------------------------------------------
# ARC-LENGTH CONTINUATION (implemented explicitly, predictor-corrector)
#
# Along the curve f(X,k)=0:  f_X*dX + f_k*dk = 0  =>  h = dX/dk = -f_k/f_X.
# Parameterize by arc length s with  ds^2 = dk^2 + dX^2  and  dX = h*dk.
# The unit tangent (t_k,t_X) satisfies f_k*t_k + f_X*t_X = 0, i.e.
#   (t_k, t_X) ~ (f_X, -f_k) = (f_X, X),  then normalized so t_k^2+t_X^2=1.
# Near a fold f_X->0 so h=dX/dk diverges; but then t_k->0 and t_X->+-1,
# so the *arc-length* step ds stays finite -> we glide through the fold.
# ----------------------------------------------------------------------
def tangent(X, k, prev):
    tk, tX = fX(X, k), X            # (f_X, -f_k)
    norm = np.hypot(tk, tX)
    tk, tX = tk / norm, tX / norm
    if tk * prev[0] + tX * prev[1] < 0:   # keep a consistent travel direction
        tk, tX = -tk, -tX
    return tk, tX

# start on the curve: pick small X, get matching k from f=0  (k = F(X)/X)
X = 3.0
k = (g0 + g1 * hill(X)) / X
prev_t = (-1.0, 1.0)                # initial heading: k decreasing, X increasing

ds, n_steps = 0.5, 1400
ks, Xs = [k], [X]
for _ in range(n_steps):
    tk, tX = tangent(X, k, prev_t)
    prev_t = (tk, tX)
    # ---- predictor: step ds along the unit tangent (finite even at folds) ----
    k0, X0 = k, X
    kp, Xp = k0 + ds * tk, X0 + ds * tX
    # ---- corrector: Newton on {f=0, pseudo-arclength constraint} ----
    kc, Xc = kp, Xp
    for _ in range(50):
        R1 = f(Xc, kc)
        R2 = (Xc - X0) * tX + (kc - k0) * tk - ds
        # Jacobian rows: [df/dX, df/dk] and [tX, tk]
        a, b = fX(Xc, kc), fk(Xc, kc)
        det = a * tk - b * tX
        dX = (R1 * tk - b * R2) / det
        dk = (a * R2 - R1 * tX) / det
        Xc, kc = Xc - dX, kc - dk
        if abs(dX) + abs(dk) < 1e-10:
            break
    X, k = Xc, kc
    ks.append(k); Xs.append(X)
    if X > 320 or k <= 0:          # traced the whole curve
        break
ks, Xs = np.array(ks), np.array(Xs)

# ----------------------------------------------------------------------
# Fold locations: dk/dX=0 on k=F(X)/X  <=>  f_X=0 on the curve.
# Scan the traced tangent's t_k for sign changes (t_k=0 at a fold).
# ----------------------------------------------------------------------
fX_vals = np.array([fX(x, kk) for x, kk in zip(Xs, ks)])
fold_idx = np.where(np.sign(fX_vals[:-1]) != np.sign(fX_vals[1:]))[0]
print("Number of arc-length points traced:", len(Xs))
print("k range covered by arc-length trace: [%.4f, %.4f]" % (ks.min(), ks.max()))
print("X range covered by arc-length trace: [%.4f, %.4f]" % (Xs.min(), Xs.max()))
print("Number of folds detected (t_k sign changes):", len(fold_idx))
for j, i in enumerate(fold_idx, 1):
    print("Fold %d near: k = %.5f , X = %.5f" % (j, ks[i], Xs[i]))
if len(fold_idx) == 2:
    klo, khi = sorted([ks[fold_idx[0]], ks[fold_idx[1]]])
    print("Bistable window in k: (%.5f, %.5f)" % (klo, khi))

# ----------------------------------------------------------------------
# CHECK: plain k-continuation (step k, solve f=0 for X by Newton from the
# previous X). Start on the upper branch and march k upward. At the upper
# fold the upper solution vanishes, so Newton cannot follow it and jumps
# to the far branch -> the unstable middle is never traced (it "stalls").
# ----------------------------------------------------------------------
kgrid = np.linspace(0.20, 0.60, 400)
Xk = np.full_like(kgrid, np.nan)
Xcur = 260.0                        # start high on the upper branch
stall_k = None
for i, kk in enumerate(kgrid):
    x = Xcur
    for _ in range(100):            # Newton in X only (k fixed)
        d = fX(x, kk)
        x -= f(x, kk) / d
    if i > 0 and abs(x - Xcur) > 40.0 and stall_k is None:
        stall_k = kgrid[i - 1]      # discontinuous branch jump = stall point
    Xk[i] = x
    Xcur = x
print("Plain k-continuation stalls (upper branch lost) near k =",
      "%.5f" % stall_k if stall_k is not None else "n/a")

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.plot(ks, Xs, '-', lw=2, color='C0',
         label='arc-length trace (full S-curve)')
plt.plot(kgrid, Xk, '--', lw=1.3, color='C3',
         label='plain k-continuation (jumps at fold)')
for j, i in enumerate(fold_idx):
    plt.plot(ks[i], Xs[i], 'ks', ms=8,
             label='fold' if j == 0 else None)
plt.xlabel('control parameter k')
plt.ylabel('equilibrium X')
plt.title('S-curve of self-activating gene via arc-length continuation')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.6.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: the check confirms success because the arc-length curve "
      "passes continuously through both folds and covers the middle "
      "(unstable) branch that plain k-continuation skips by jumping, so "
      "recovering that connected S-shape proves the fold singularity was "
      "handled rather than stalled at.")
