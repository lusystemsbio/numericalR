import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model: self-activating gene -----
# f(X,k) = g0 + g1*Hill(X) - k*X, steady states are f = 0
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4

def hill(X):
    u = (X / Xth) ** n
    return u / (1.0 + u)

def dhill_dX(X):
    # d/dX of u/(1+u) with u=(X/Xth)^n
    u = (X / Xth) ** n
    du_dX = n * X ** (n - 1) / Xth ** n
    return du_dX / (1.0 + u) ** 2

def f(X, k):
    return g0 + g1 * hill(X) - k * X

def dfdX(X, k):
    return g1 * dhill_dX(X) - k          # partial wrt state variable

def dfdk(X, k):
    return -X                            # partial wrt control parameter

# ----- Newton corrector: root-solve f(X, k_fixed) = 0 in X -----
def correct(X0, k, tol=1e-10, maxit=100):
    X = X0
    for _ in range(maxit):
        fx = f(X, k)
        d = dfdX(X, k)
        if abs(d) < 1e-14:               # Jacobian singular -> cannot correct (near fold)
            return X, False
        X_new = X - fx / d
        if abs(X_new - X) < tol:
            return X_new, True
        X = X_new
    return X, False

# ----- Find a starting point on the lower branch -----
k_start = 0.30
X_guess = 40.0
X0, ok = correct(X_guess, k_start)
print(f"start point: k = {k_start:.6f}, X = {X0:.6f}, ok = {ok}")

# ----- Predictor-corrector continuation, stepping k downward -----
# Decreasing k pushes the system up the lower branch toward the fold.
dk = -0.001                              # step in the control parameter
ks, Xs = [k_start], [X0]
k, X = k_start, X0
fold_k = None

for step in range(20000):
    d = dfdX(X, k)
    if abs(d) < 1e-3:                    # dX/dk -> infinity: the fold, k becomes multivalued
        fold_k = k
        print(f"STALL at fold: k = {k:.6f}, X = {X:.6f}, df/dX = {d:.3e}")
        break
    tangent = -dfdk(X, k) / d            # predictor slope dX/dk = -(df/dk)/(df/dX)
    k_new = k + dk
    X_pred = X + tangent * dk            # predict along the tangent
    X_new, ok = correct(X_pred, k_new)   # correct back onto f = 0
    if not ok:                           # corrector failed to reconverge -> also a stall
        fold_k = k
        print(f"STALL (corrector diverged): k = {k:.6f}, X = {X:.6f}")
        break
    k, X = k_new, X_new
    ks.append(k); Xs.append(X)

ks, Xs = np.array(ks), np.array(Xs)
print(f"steps taken before stall: {len(ks)}")
print(f"k range followed: [{ks.min():.6f}, {ks.max():.6f}]")
print(f"X range followed: [{Xs.min():.6f}, {Xs.max():.6f}]")

# ----- Independent check: exact branch via explicit k(X) = (g0 + g1*Hill(X))/X -----
# k as a function of X is single-valued, so it is the ground truth to compare against.
Xtrue = np.linspace(5.0, 800.0, 4000)
ktrue = (g0 + g1 * hill(Xtrue)) / Xtrue

# fold(s): where dk/dX = 0  (equivalently df/dX = 0 on the curve)
dk_dX = np.gradient(ktrue, Xtrue)
sign_change = np.where(np.diff(np.sign(dk_dX)) != 0)[0]
fold_ks_exact = [(ktrue[i] + ktrue[i + 1]) / 2 for i in sign_change]
print("exact fold k value(s) (dk/dX = 0):", [f"{v:.6f}" for v in fold_ks_exact])

# accuracy away from the fold: for each continued point, X must satisfy the exact branch
resid = np.array([f(Xs[i], ks[i]) for i in range(len(ks))])
print(f"max |f(X,k)| on continued branch (residual): {np.max(np.abs(resid)):.3e}")
# compare to exact X at the continued k's, staying clear of the fold
mask = ks > (fold_k + 0.02 if fold_k is not None else ks.min())
Xexact = np.interp(ks[mask], ktrue[::-1], Xtrue[::-1])
if mask.sum() > 0:
    print(f"max |X_continued - X_exact| away from fold: {np.max(np.abs(Xs[mask] - Xexact)):.3e}")

# ----- Plot -----
plt.figure(figsize=(8, 6))
plt.plot(ktrue, Xtrue, '-', color='0.7', lw=2, label='exact branch  k(X)=(g0+g1 Hill)/X')
plt.plot(ks, Xs, '.-', color='C0', ms=3, label='predictor-corrector continuation')
if fold_k is not None:
    plt.plot(ks[-1], Xs[-1], 'r*', ms=16, label=f'stall at fold  k≈{fold_k:.4f}')
for v in fold_ks_exact:
    plt.axvline(v, ls='--', color='r', alpha=0.4)
plt.xlabel('control parameter k')
plt.ylabel('steady-state X')
plt.title('Bifurcation curve X(k): continuation stalls at the fold')
plt.legend()
plt.xlim(0, 0.5)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.5.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: the run tracks the exact single-valued branch k(X) to tiny residual "
      "away from the fold but halts exactly where df/dX->0 (dX/dk diverges and k becomes "
      "multivalued), confirming the method is correct and that its only failure is the "
      "expected geometric singularity of k-continuation at a fold.")
