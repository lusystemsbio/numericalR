import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Model: self-activating gene
#   f(X,k) = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X
# Equilibria are the roots f(X,k)=0. As a function of k this set is
# S-shaped (two stable branches + an unstable middle) with two folds
# where the tangent goes vertical (dX/dk -> infinity). Plain
# k-continuation stalls there; arc-length continuation does not.
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def hill(X):
    u = (X / Xth)**n
    return u / (1.0 + u)

def hill_prime(X):
    # dH/dX with u=(X/Xth)^n, du/dX = n*u/X
    u = (X / Xth)**n
    return (n * u / X) / (1.0 + u)**2

def f(X, k):
    return g0 + g1 * hill(X) - k * X

def fX(X, k):
    return g1 * hill_prime(X) - k          # partial df/dX

def fk(X, k):
    return -X                              # partial df/dk

# ---------------------------------------------------------------
# Arc-length continuation, done explicitly.
# We parameterize the curve by arc length s, with a state (X,k).
# A step of length ds is split into (dX, dk) obeying
#     ds^2 = dk^2 + dX^2,   dX = h*dk   (h = dX/dk = tangent slope).
# Because h can diverge at a fold, we never invert it: instead we
# work with the *unit tangent* (tX,tk) = (dX/ds, dk/ds), which stays
# finite (tk -> 0 at a fold, tX -> +/-1). This keeps every step finite.
# ---------------------------------------------------------------

# Start on the lower branch at a known equilibrium: pick X0, get k0 exactly.
X, k = 20.0, 0.0
k = (g0 + g1 * hill(X)) / X        # exact root of f(X,k)=0 at this X

# Initial unit tangent: differentiate f=0 -> fX*tX + fk*tk = 0,
# plus tX^2+tk^2=1. Choose the sign that moves toward larger X.
def unit_tangent(X, k, prev=None):
    a, b = fX(X, k), fk(X, k)               # fX*tX + fk*tk = 0
    tX, tk = -b, a                          # a vector orthogonal to (a,b)
    norm = np.hypot(tX, tk)
    tX, tk = tX / norm, tk / norm
    if prev is None:
        if tX < 0:                          # start heading to larger X
            tX, tk = -tX, -tk
    else:                                   # keep direction consistent
        if tX * prev[0] + tk * prev[1] < 0:
            tX, tk = -tX, -tk
    return tX, tk

ds = 2.0
Xs, ks, tks = [X], [k], []
tX, tk = unit_tangent(X, k)
tks.append(tk)

for _ in range(4000):
    # ---- Predictor: linear step along the tangent (finite even at folds)
    X0, k0 = X, k
    Xp, kp = X0 + ds * tX, k0 + ds * tk

    # ---- Corrector: Newton on the augmented system
    #   f(X,k) = 0
    #   (X-X0)*tX + (k-k0)*tk - ds = 0   (pseudo-arclength constraint)
    Xc, kc = Xp, kp
    for _ in range(50):
        F1 = f(Xc, kc)
        F2 = (Xc - X0) * tX + (kc - k0) * tk - ds
        J = np.array([[fX(Xc, kc), fk(Xc, kc)],
                      [tX,          tk       ]])
        dX_, dk_ = np.linalg.solve(J, [-F1, -F2])
        Xc += dX_
        kc += dk_
        if abs(dX_) + abs(dk_) < 1e-10:
            break

    X, k = Xc, kc
    Xs.append(X); ks.append(k)

    # ---- New tangent for next step, sign kept consistent
    tX, tk = unit_tangent(X, k, prev=(tX, tk))
    tks.append(tk)

    # Stop once we have climbed well onto the upper branch
    if X > 3000:
        break

Xs, ks, tks = np.array(Xs), np.array(ks), np.array(tks)

# ---------------------------------------------------------------
# Locate the folds: a fold is where the k-tangent tk = dk/ds changes
# sign (curve momentarily vertical, dX/dk -> infinity).
# ---------------------------------------------------------------
fold_idx = np.where(np.sign(tks[:-1]) != np.sign(tks[1:]))[0]
fold_ks = [0.5 * (ks[i] + ks[i + 1]) for i in fold_idx]
fold_Xs = [0.5 * (Xs[i] + Xs[i + 1]) for i in fold_idx]

print(f"Number of traced points: {len(Xs)}")
print(f"Number of folds detected (dk/ds sign changes): {len(fold_idx)}")
for j, (kf, Xf) in enumerate(zip(fold_ks, fold_Xs), 1):
    print(f"Fold {j}: k = {kf:.6f}, X = {Xf:.6f}")

k_low, k_high = min(fold_ks), max(fold_ks)
print(f"Bistable k-range (between folds): k_low = {k_low:.6f}, k_high = {k_high:.6f}")

# ---------------------------------------------------------------
# Check: at a control value inside the bistable window the curve must
# cross that k three times (lower stable, unstable middle, upper stable).
# Count crossings of the traced curve with k = k_mid.
# ---------------------------------------------------------------
k_mid = 0.5 * (k_low + k_high)
d = ks - k_mid
cross = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
X_at_kmid = []
for i in cross:
    # linear interpolation in X across the crossing
    t = (k_mid - ks[i]) / (ks[i + 1] - ks[i])
    X_at_kmid.append(Xs[i] + t * (Xs[i + 1] - Xs[i]))
X_at_kmid.sort()

print(f"Sample control value k_mid = {k_mid:.6f}")
print(f"Number of equilibria found at k_mid: {len(X_at_kmid)}")
for j, Xv in enumerate(X_at_kmid, 1):
    print(f"  equilibrium {j}: X = {Xv:.6f}")

# Stability of each equilibrium at k_mid: stable if fX < 0, unstable if fX > 0
for j, Xv in enumerate(X_at_kmid, 1):
    slope = fX(Xv, k_mid)
    print(f"  equilibrium {j} fX = {slope:.6e} -> {'stable' if slope < 0 else 'unstable'}")

smooth = (len(fold_idx) == 2) and (len(X_at_kmid) == 3)
print(f"Traced smoothly through both folds and recovered full S-shape: {smooth}")

# ---------------------------------------------------------------
# Plot the full S-curve
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(ks, Xs, '-', color='steelblue', lw=1.8, label='arc-length continuation')
plt.plot(fold_ks, fold_Xs, 'ko', ms=8, label='folds')
plt.axvline(k_mid, color='gray', ls='--', lw=1, label=f'k = {k_mid:.4f}')
plt.plot([k_mid] * len(X_at_kmid), X_at_kmid, 'rs', ms=7,
         label=f'{len(X_at_kmid)} equilibria at k')
plt.xlabel('control parameter k')
plt.ylabel('equilibrium X')
plt.title('Self-activating gene: full S-curve via arc-length continuation')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.6.1_s1.png")

# One-sentence explanation of why the check confirms the result:
# Finding exactly two sign-changes of dk/ds (the two folds) plus three
# equilibria at a single k inside the fold window proves the tracer walked
# continuously from one stable branch, through both vertical folds and the
# unstable middle, onto the other stable branch -- precisely the traversal
# where plain k-continuation stalls because dX/dk diverges.
print("Check rationale: two dk/ds sign-changes and three coexisting equilibria at one k "
      "prove continuous passage through both vertical folds and the unstable middle branch, "
      "the exact traversal that plain k-continuation cannot make because dX/dk diverges there.")
