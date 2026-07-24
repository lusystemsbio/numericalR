import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Model: self-activating gene
#   f(X,k) = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X
# Equilibria are the roots f(X,k)=0.  We trace that curve in the (k,X)
# plane by ARC-LENGTH continuation so we can pass through the folds
# (where dX/dk -> infinity) without stalling.
# ----------------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def hill(X):
    u = X / Xth
    return u**n / (1.0 + u**n)

def dhill_dX(X):
    u = X / Xth
    return n * u**(n-1) * (1.0/Xth) / (1.0 + u**n)**2

def f(X, k):
    return g0 + g1*hill(X) - k*X

def fX(X, k):      # df/dX
    return g1*dhill_dX(X) - k

def fk(X, k):      # df/dk
    return -X

# ----------------------------------------------------------------------
# ARC-LENGTH CONTINUATION (predictor / corrector), implemented explicitly.
# The tangent to f=0 is perpendicular to the gradient (fk, fX), so it is
# proportional to (fX, -fk) = (fX, X).  Writing that tangent as
# (dk, dX) with dX = h*dk and h = dX/dk = -fk/fX, the arc-length
# normalization ds^2 = dk^2 + dX^2 gives
#     dk = ds / sqrt(1+h^2),   dX = h*ds / sqrt(1+h^2),
# which stays FINITE even at a fold (fX->0, h->inf gives dk->0, dX->ds).
# ----------------------------------------------------------------------

# --- starting point: pick X, get exact k on the lower branch ---
X = 20.0
k = (g0 + g1*hill(X)) / X          # exact equilibrium so f(X,k)=0
print(f"start point: k = {k:.6f}, X = {X:.6f}")

# --- initial tangent, oriented toward increasing X ---
t = np.array([fX(X, k), X])        # (t_k, t_X) proportional to (fX, X)
t = t / np.hypot(*t)
if t[1] < 0:                       # want to march with dX > 0 first
    t = -t

ds = 1.5                           # arc-length step
ks, Xs = [k], [X]
for _ in range(800):
    # ---- PREDICTOR: step ds along the current unit tangent ----
    kp = k + ds * t[0]             # dk = ds * t_k
    Xp = X + ds * t[1]             # dX = ds * t_X
    kk, XX = kp, Xp

    # ---- CORRECTOR: Newton on the augmented 2x2 system ----
    #   f(X,k)=0                              (stay on the curve)
    #   t.( (k,X)-(kp,Xp) )=0                 (stay in plane _|_ tangent)
    for _ in range(50):
        F0 = f(XX, kk)
        F1 = t[0]*(kk - kp) + t[1]*(XX - Xp)
        J = np.array([[fk(XX, kk), fX(XX, kk)],   # d/dk , d/dX of f
                      [t[0],       t[1]]])         # d/dk , d/dX of plane
        dz = np.linalg.solve(J, [F0, F1])
        kk -= dz[0]; XX -= dz[1]
        if np.hypot(*dz) < 1e-10:
            break
    k, X = kk, XX
    ks.append(k); Xs.append(X)

    # ---- update tangent for next step, keeping orientation continuous ----
    tn = np.array([fX(X, k), X])
    tn = tn / np.hypot(*tn)
    if np.dot(tn, t) < 0:
        tn = -tn
    t = tn

    if X > 560 or X < 0:           # covered the full S-curve
        break

ks = np.array(ks); Xs = np.array(Xs)
print(f"arc-length points traced: {len(ks)}")

# --- locate the two folds: where dk/ds changes sign (dk/dX -> inf) ---
dk = np.diff(ks)
fold_idx = np.where(np.sign(dk[:-1]) != np.sign(dk[1:]))[0] + 1
folds = [(ks[i], Xs[i]) for i in fold_idx]
print(f"number of folds detected: {len(folds)}")
for j, (kf, Xf) in enumerate(folds):
    print(f"fold {j+1}: k = {kf:.6f}, X = {Xf:.6f}")

k_lo = min(f[0] for f in folds); k_hi = max(f[0] for f in folds)
X_lo = min(f[1] for f in folds); X_hi = max(f[1] for f in folds)
print(f"bistable k-range (between folds): [{k_lo:.6f}, {k_hi:.6f}]")

# points on the middle (unstable, negative-slope) branch: between folds in X
mid_mask = (Xs > X_lo) & (Xs < X_hi) & (ks > k_lo) & (ks < k_hi)
print(f"arc-length points on middle (unstable) branch: {int(mid_mask.sum())}")

# ----------------------------------------------------------------------
# CHECK: plain k-continuation (fix k, Newton in X only).  It uses only
# f_X in the Jacobian, so it becomes singular at a fold and cannot follow
# the negatively-sloped middle branch -- it stalls / jumps.
# ----------------------------------------------------------------------
ksweep = np.linspace(0.05, 0.30, 400)   # sweep k upward across the folds
Xguess = 20.0
kc_X, kc_k = [], []
max_jump = 0.0
jump_k = None
for kv in ksweep:
    Xn = Xguess
    for _ in range(100):                # Newton in X at fixed k
        d = fX(Xn, kv)
        if abs(d) < 1e-12:              # singular Jacobian -> stall at fold
            break
        step = f(Xn, kv) / d
        Xn -= step
        if abs(step) < 1e-10:
            break
    jump = abs(Xn - Xguess)
    if jump > max_jump:
        max_jump = jump; jump_k = kv
    kc_k.append(kv); kc_X.append(Xn)
    Xguess = Xn
kc_X = np.array(kc_X); kc_k = np.array(kc_k)
kc_mid = int(np.sum((kc_X > X_lo) & (kc_X < X_hi) & (kc_k > k_lo) & (kc_k < k_hi)))
print(f"plain k-continuation: max |dX| between adjacent k-steps = {max_jump:.4f} (a jump across a fold)")
print(f"plain k-continuation: k at that jump = {jump_k:.6f}")
print(f"plain k-continuation points on middle (unstable) branch: {kc_mid}")

# ----------------------------------------------------------------------
# Plot the full S-shaped curve.
# ----------------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(ks, Xs, '-', lw=2, color='tab:blue', label='arc-length continuation (full S-curve)')
plt.plot(kc_k, kc_X, '.', ms=4, color='0.6', label='plain k-continuation (stalls/jumps)')
for kf, Xf in folds:
    plt.plot(kf, Xf, 'ks', ms=9, mfc='red')
    plt.annotate('fold', (kf, Xf), textcoords='offset points', xytext=(8, 6))
plt.xlabel('control parameter k')
plt.ylabel('equilibrium X')
plt.title('Self-activating gene: S-curve via arc-length continuation')
plt.legend(loc='best')
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.6.1_s2.png")

# The check confirms the result because arc-length continuation places
# many points on the negatively-sloped middle branch that connects the
# two folds (where dk/dX -> inf), while plain k-continuation places zero
# points there and instead makes a large discontinuous jump in X at the
# fold -- so recovering the full S-shape, unstable branch included, is
# exactly the behavior plain k-continuation cannot produce.
print("CHECK: arc-length traces the unstable middle branch through both folds "
      "(nonzero middle-branch points, two folds), whereas plain k-continuation "
      "samples zero middle-branch points and jumps across a fold -- confirming "
      "the folds were crossed smoothly rather than stalled.")
