import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Model: self-activating gene at steady state f(X,k)=0
#   f = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def hill(X):
    u = (X / Xth)**n
    return u / (1.0 + u)

def dhill(X):                       # d/dX of Hill term
    u = (X / Xth)
    # derivative of u^n/(1+u^n) w.r.t X, chain rule with u=X/Xth
    return n * u**(n - 1) * (1.0 / Xth) / (1.0 + u**n)**2

def f(X, k):                        # residual
    return g0 + g1 * hill(X) - k * X

def fX(X, k):                       # partial df/dX  (its sign = stability)
    return g1 * dhill(X) - k

def fk(X, k):                       # partial df/dk
    return -X

# ---------------------------------------------------------------
# Arc-length (pseudo-arclength) continuation, implemented by hand.
# We treat BOTH k and X as unknowns and advance by arc length s.
# Each step splits into (dk, dX) constrained by  ds^2 = dk^2 + dX^2.
# Along the curve f=0:  fX*dX + fk*dk = 0  =>  dX/dk = -fk/fX = h.
# At a fold fX->0 so h=dX/dk diverges; stepping in k alone stalls,
# but the unit tangent (dk,dX) stays finite, so arc length does not.
# ---------------------------------------------------------------

# --- starting point on the lower branch: pick X, solve f=0 for k ---
X = 1.0
k = (g0 + g1 * hill(X)) / X         # exact k making f(X,k)=0 at this X
k0curve = k

# --- initial unit tangent, oriented toward increasing X ---
# tangent is perpendicular to gradient (fk,fX): direction (fX, -fk)=(fX, X)
tk, tX = fX(X, k), -fk(X, k)
nrm = np.hypot(tk, tX)
tk, tX = tk / nrm, tX / nrm
if tX < 0:                          # force initial motion to raise X
    tk, tX = -tk, -tX

ds = 3.0                            # arc-length step
ks, Xs, stab = [k], [X], [fX(X, k) < 0]

for step in range(600):
    k0, X0 = k, X
    tk0, tX0 = tk, tX

    # ---- predictor: move one arc-length along the tangent ----
    kp, Xp = k0 + ds * tk0, X0 + ds * tX0

    # ---- corrector: Newton on augmented 2x2 system ----
    #   F1 = f(X,k) = 0
    #   F2 = tk*(k-k0) + tX*(X-X0) - ds = 0   (arc-length constraint)
    k, X = kp, Xp
    for _ in range(50):
        F1 = f(X, k)
        F2 = tk0 * (k - k0) + tX0 * (X - X0) - ds
        # Jacobian rows: d/dk, d/dX
        J = np.array([[fk(X, k), fX(X, k)],
                      [tk0,      tX0     ]])
        dk, dX = np.linalg.solve(J, [-F1, -F2])
        k += dk
        X += dX
        if abs(dk) + abs(dX) < 1e-10:
            break

    # ---- update tangent, keep it pointing the same way (no U-turn) ----
    tk, tX = fX(X, k), -fk(X, k)
    nrm = np.hypot(tk, tX)
    tk, tX = tk / nrm, tX / nrm
    if tk * tk0 + tX * tX0 < 0:
        tk, tX = -tk, -tX

    ks.append(k)
    Xs.append(X)
    stab.append(fX(X, k) < 0)       # stable where df/dX < 0

    if X > 1200 or X < 0:           # traced past the upper branch
        break

ks, Xs, stab = np.array(ks), np.array(Xs), np.array(stab)

# ---------------------------------------------------------------
# Fold detection: a fold is where the tangent's k-component changes
# sign (dk/ds flips), i.e. k reverses direction while X keeps going.
# ---------------------------------------------------------------
dk_ds = np.gradient(ks)             # sign of progress in k along the curve
sign_changes = np.where(np.diff(np.sign(dk_ds)) != 0)[0]
folds = [(ks[i], Xs[i]) for i in sign_changes]

# ---------------------------------------------------------------
# Bistability check: for a k inside the fold interval, count how many
# points of the traced curve lie near it (should be 3 -> full S).
# ---------------------------------------------------------------
if len(folds) >= 2:
    kf = sorted([folds[0][0], folds[1][0]])
    k_lo, k_hi = kf[0], kf[1]
    k_test = 0.5 * (k_lo + k_hi)
    # crossings of the curve with the vertical line k = k_test
    crossings = np.where(np.diff(np.sign(ks - k_test)) != 0)[0]
    n_branches = len(crossings)
else:
    k_lo = k_hi = k_test = float('nan')
    n_branches = 0

# ---------------------------------------------------------------
# Print numerical results
# ---------------------------------------------------------------
print(f"Model params: g0={g0}, g1={g1}, Xth={Xth}, n={n}")
print(f"Number of continuation points traced: {len(ks)}")
print(f"Starting point (k, X): ({k0curve:.6f}, 1.000000)")
print(f"Final point (k, X): ({ks[-1]:.6f}, {Xs[-1]:.6f})")
print(f"Number of folds detected: {len(folds)}")
for j, (kf_, Xf_) in enumerate(folds, 1):
    print(f"Fold {j}: k = {kf_:.6f}, X = {Xf_:.6f}")
print(f"Bistable k-interval: [{k_lo:.6f}, {k_hi:.6f}]")
print(f"Bistable interval width in k: {k_hi - k_lo:.6f}")
print(f"Test k inside interval: {k_test:.6f}")
print(f"Number of coexisting X-branches at test k: {n_branches}")
print(f"Stable points (df/dX<0): {int(np.sum(stab))}")
print(f"Unstable points (df/dX>0): {int(np.sum(~stab))}")

# ---------------------------------------------------------------
# Plot the full S-curve
# ---------------------------------------------------------------
plt.figure(figsize=(8, 6))
stable_mask = stab
plt.plot(ks[stable_mask], Xs[stable_mask], '.', color='C0', ms=4, label='stable (df/dX<0)')
plt.plot(ks[~stable_mask], Xs[~stable_mask], '.', color='C3', ms=4, label='unstable (df/dX>0)')
plt.plot(ks, Xs, '-', color='gray', lw=0.6, alpha=0.6, zorder=0)
for j, (kf_, Xf_) in enumerate(folds, 1):
    plt.plot(kf_, Xf_, 'ks', ms=8, mfc='yellow', label='fold' if j == 1 else None)
if n_branches:
    plt.axvline(k_test, color='green', ls='--', lw=1, label=f'k={k_test:.3f} ({n_branches} states)')
plt.xlabel('control parameter k')
plt.ylabel('steady-state X')
plt.title('Self-activating gene: full S-curve via arc-length continuation')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.6.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Check explanation: detecting two folds where dk/ds reverses sign while X "
      "advances monotonically, together with three coexisting X-values at a single "
      "k inside the fold interval, confirms the arc-length method rounded both turning "
      "points and recovered the full S (both stable branches plus the unstable middle) "
      "exactly where k-only continuation would reverse and stall.")
