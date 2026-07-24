import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Model: self-activating gene
#   f(X,k) = g0 + g1*(X/Xth)^n / (1 + (X/Xth)^n) - k*X
# Fixed points solve f(X,k)=0. As k varies this set is an S-shaped curve
# in the (k,X) plane with two folds (turning points) => bistability.
# ----------------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(X, k):
    u = (X / Xth) ** n
    return g0 + g1 * u / (1.0 + u) - k * X

def fX(X, k):
    # partial derivative wrt X (Jacobian entry that vanishes AT the folds)
    u = (X / Xth) ** n
    return g1 * n * u / (X * (1.0 + u) ** 2) - k

def fk(X, k):
    # partial derivative wrt k
    return -X

# ----------------------------------------------------------------------
# Reference: because k can be written explicitly as k(X)=(g0+g1*H(X))/X,
# the curve is single valued in X, and its folds are the extrema of k(X).
# We use this only to LOCATE the folds for reporting/checking.
# ----------------------------------------------------------------------
Xg = np.linspace(1.0, 2000.0, 400000)
ug = (Xg / Xth) ** n
kg = (g0 + g1 * ug / (1.0 + ug)) / Xg
dk = np.diff(kg)
turn_idx = np.where(np.diff(np.sign(dk)) != 0)[0] + 1  # local extrema of k(X)
folds = [(kg[i], Xg[i]) for i in turn_idx]
folds.sort()  # by k
k_lo_fold, X_at_lo = folds[0]
k_hi_fold, X_at_hi = folds[-1]
print(f"Lower fold (turning point):  k = {k_lo_fold:.6f}, X = {X_at_lo:.4f}")
print(f"Upper fold (turning point):  k = {k_hi_fold:.6f}, X = {X_at_hi:.4f}")
print(f"Bistable k-range (three coexisting fixed points): "
      f"[{k_lo_fold:.6f}, {k_hi_fold:.6f}]")

# ----------------------------------------------------------------------
# Starting point on the lower (low-expression) branch: pick k, solve for X.
# ----------------------------------------------------------------------
def solve_X_at_k(k, X0):
    X = X0
    for _ in range(200):
        step = f(X, k) / fX(X, k)
        X -= step
        if abs(step) < 1e-12:
            break
    return X

k_start = 0.5
X_start = solve_X_at_k(k_start, 15.0)
print(f"Start point on lower branch: k = {k_start:.6f}, X = {X_start:.6f}")

# ======================================================================
# ARC-LENGTH (pseudo-arclength) CONTINUATION
# Idea: parameterize the curve by arc length s, not by k. A step of
# length ds is split into a (dk, dX) pair with  ds^2 = dk^2 + dX^2 and
# dX = h*dk where h = dX/dk = -fk/fX. Away from folds this is ordinary
# k-stepping; where dX/dk -> infinity (fk/fX -> inf, i.e. fX -> 0) the
# split keeps dk -> 0 while dX stays finite, so the step never blows up.
# Equivalently the unit tangent (tk,tX) is proportional to (fX, -fk),
# since fX*tX + fk*tk = 0 must hold along f=0.
# ======================================================================
ds = 2.0                 # fixed arc-length step
X, k = X_start, k_start
tk_prev, tX_prev = 0.0, 1.0   # want X increasing along the curve initially

ks, Xs = [k], [X]
max_steps = 5000
for _ in range(max_steps):
    # --- tangent from  fX*tX + fk*tk = 0  =>  (tk,tX) ~ (fX, -fk) -----
    a, b = fX(X, k), fk(X, k)          # a=fX, b=fk
    tk, tX = a, -b                     # tangent direction (unnormalized)
    norm = np.hypot(tk, tX)
    tk, tX = tk / norm, tX / norm
    # keep marching the same way along the curve (sign continuity)
    if tk * tk_prev + tX * tX_prev < 0.0:
        tk, tX = -tk, -tX
    tk_prev, tX_prev = tk, tX

    # note: h = dX/dk = tX/tk diverges at a fold, but the split below
    # uses ds directly, so dk = ds*tk -> 0 there instead of exploding.
    kp = k + ds * tk                   # predictor:  dk = ds*tk
    Xp = X + ds * tX                   #             dX = ds*tX = h*dk

    # --- Newton corrector back onto f=0 with arc-length constraint ----
    #   f(Xc,kc) = 0
    #   tX*(Xc-X) + tk*(kc-k) - ds = 0   (stay ds away along the tangent)
    Xc, kc = Xp, kp
    for _ in range(100):
        r1 = f(Xc, kc)
        r2 = tX * (Xc - X) + tk * (kc - k) - ds
        # 2x2 Jacobian: rows d(r1),d(r2) wrt (Xc,kc)
        J11, J12 = fX(Xc, kc), fk(Xc, kc)
        J21, J22 = tX, tk
        det = J11 * J22 - J12 * J21     # stays ~ X != 0 even where fX->0
        dX_ = (r1 * J22 - r2 * J12) / det
        dk_ = (J11 * r2 - J21 * r1) / det
        Xc -= dX_
        kc -= dk_
        if abs(dX_) + abs(dk_) < 1e-11:
            break

    X, k = Xc, kc
    ks.append(k)
    Xs.append(X)
    # full traverse: stop once we are far along the upper branch
    if X > 900.0 or k < 1e-3:
        break

ks = np.array(ks)
Xs = np.array(Xs)
print(f"Arc-length continuation: traced {len(Xs)} points")
print(f"  k spanned: [{ks.min():.6f}, {ks.max():.6f}]")
print(f"  X spanned: [{Xs.min():.4f}, {Xs.max():.4f}]")
print(f"  End point: k = {ks[-1]:.6f}, X = {Xs[-1]:.4f}")

# ======================================================================
# CHECK: plain k-continuation (step k, Newton-solve X each time).
# It must stall at the lower fold, where fX -> 0 makes the Newton
# derivative singular and no nearby same-branch solution exists.
# ======================================================================
kc = k_start
Xc = X_start
k_seq = []
X_seq = []
dk_step = -0.002          # march k downward off the lower branch
stall_k = None
prevX = Xc
for _ in range(100000):
    kc += dk_step
    if kc <= 0:
        break
    # Newton in X only (this is the "one-variable" k-continuation)
    ok = True
    Xn = Xc
    for _ in range(50):
        d = fX(Xn, kc)
        if abs(d) < 1e-6:          # Jacobian singular near the fold
            ok = False
            break
        s = f(Xn, kc) / d
        Xn -= s
        if abs(s) < 1e-10:
            break
    if (not ok) or (not np.isfinite(Xn)) or abs(Xn - prevX) > 50.0:
        stall_k = kc - dk_step     # last k that still worked
        break
    Xc, prevX = Xn, Xn
    k_seq.append(kc)
    X_seq.append(Xc)

if stall_k is not None:
    print(f"Plain k-continuation STALLED at k = {stall_k:.6f} "
          f"(near lower fold k = {k_lo_fold:.6f})")
    print(f"  it only reached X = {X_seq[-1] if X_seq else X_start:.4f}, "
          f"covering a single branch")
else:
    print("Plain k-continuation did not stall (unexpected)")

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.plot(ks, Xs, '-', color='tab:blue', lw=2,
         label='arc-length continuation (full S-curve)')
if X_seq:
    plt.plot(k_seq, X_seq, '--', color='tab:red', lw=2,
             label='plain k-continuation (stalls)')
plt.plot([k_lo_fold, k_hi_fold], [X_at_lo, X_at_hi], 'ks',
         ms=8, label='folds')
if stall_k is not None:
    plt.axvline(stall_k, color='tab:red', ls=':', alpha=0.6)
plt.xlabel('k (control parameter)')
plt.ylabel('X (steady-state expression)')
plt.title('Self-activating gene: S-curve via arc-length continuation')
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.6.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("Check rationale: arc-length continuation reaches X-values on both "
      "the low and high stable branches AND the unstable middle branch, "
      "sweeping k back and forth past both fold k-values, whereas plain "
      "k-continuation halts at the lower fold where fX->0; recovering the "
      "full S where the naive method stalls confirms the folds were "
      "traversed without stalling.")
