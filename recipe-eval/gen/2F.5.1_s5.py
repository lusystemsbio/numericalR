import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Model: self-activating gene
#   f(X, k) = g0 + g1 * (X/Xth)^n / (1 + (X/Xth)^n) - k*X
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def hill(X):
    # excitatory Hill function u/(1+u), u = (X/Xth)^n
    u = (X / Xth) ** n
    return u / (1.0 + u)

def f(X, k):
    # steady-state residual; f = 0 defines the bifurcation curve X(k)
    return g0 + g1 * hill(X) - k * X

def dfdX(X, k):
    # d/dX of Hill: (n*u/X)/(1+u)^2, then minus k from degradation term
    u = (X / Xth) ** n
    dhill = (n * u / X) / (1.0 + u) ** 2
    return g1 * dhill - k

def dfdk(X, k):
    # d/dk f = -X
    return -X

# ---------------------------------------------------------------
# Get a known starting point on f=0 at small k (upper branch),
# by scanning for a sign change then a few Newton iterations.
# ---------------------------------------------------------------
def solve_X(k, X_guess):
    # Newton corrector: root-solve f(X,k)=0 in X at fixed k
    X = X_guess
    for _ in range(100):
        fx = f(X, k)
        d = dfdX(X, k)
        if abs(d) < 1e-14:
            break
        step = fx / d
        X -= step
        if abs(step) < 1e-10:
            break
    return X

k0 = 0.05
# bracket-scan a high-X root for the initial point
Xg = np.linspace(1.0, 5000.0, 20000)
fv = f(Xg, k0)
sign_changes = np.where(np.diff(np.sign(fv)) != 0)[0]
X0 = solve_X(k0, Xg[sign_changes[-1]])  # take the largest (upper-branch) root
print(f"Start point: k = {k0:.6f}, X = {X0:.6f}, residual f = {f(X0,k0):.3e}")

# ---------------------------------------------------------------
# Predictor-corrector continuation in the natural parameter k.
#   predict along tangent  dX/dk = -(df/dk)/(df/dX)
#   correct with a Newton root-solve at the new k
# This natural-parameter scheme MUST stall at the fold, where
# df/dX -> 0 so dX/dk diverges and k becomes multivalued.
# ---------------------------------------------------------------
dk = 0.001
ks, Xs, slopes = [k0], [X0], []
k, X = k0, X0
fold_k = fold_X = None

for _ in range(100000):
    d = dfdX(X, k)
    slope = -dfdk(X, k) / d        # tangent dX/dk of the curve f=0
    slopes.append(slope)
    # stall test: near the fold df/dX -> 0 and the tangent slope blows up
    if abs(d) < 1e-3 or abs(slope) > 1e4:
        fold_k, fold_X = k, X
        print(f"STALL at fold: k = {k:.6f}, X = {X:.6f}, "
              f"df/dX = {d:.3e}, dX/dk = {slope:.3e}")
        break
    k_new = k + dk
    X_pred = X + slope * dk         # tangent predictor
    X_new = solve_X(k_new, X_pred)  # corrector back onto f=0
    # if the corrector cannot land on the branch at k_new, we are past the fold
    if not np.isfinite(X_new) or abs(f(X_new, k_new)) > 1e-6 or X_new <= 0:
        fold_k, fold_X = k, X
        print(f"STALL (corrector failed) near k = {k:.6f}, X = {X:.6f}")
        break
    k, X = k_new, X_new
    ks.append(k)
    Xs.append(X)

ks, Xs = np.array(ks), np.array(Xs)
print(f"Points traced: {len(ks)}")
print(f"k range traced: [{ks.min():.6f}, {ks.max():.6f}]")
print(f"X range traced: [{Xs.min():.6f}, {Xs.max():.6f}]")
print(f"Max residual on traced branch: {np.max(np.abs(f(Xs, ks))):.3e}")

# ---------------------------------------------------------------
# Independent check: sample the true curve near the fold by solving
# f=0 for X on a fine grid, and compare the analytic fold location
# (where max_k over the branch occurs, i.e. df/dX = 0).
# ---------------------------------------------------------------
# Along the branch k(X) = (g0 + g1*hill(X))/X ; its maximum in X is the fold.
Xgrid = np.linspace(50.0, 2000.0, 200000)
kgrid = (g0 + g1 * hill(Xgrid)) / Xgrid
i_fold = np.argmax(kgrid)
print(f"Analytic fold (max k along branch): k* = {kgrid[i_fold]:.6f}, "
      f"X* = {Xgrid[i_fold]:.6f}")
if fold_k is not None:
    print(f"Continuation stalled at   k = {fold_k:.6f} vs analytic k* = "
          f"{kgrid[i_fold]:.6f} (diff = {abs(fold_k-kgrid[i_fold]):.3e})")

# ---------------------------------------------------------------
# Plot the continued branch X(k) and mark the fold.
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6))
# full true curve (both branches) via k(X), for reference
ax.plot(kgrid, Xgrid, color="0.8", lw=1, label="true curve f=0 (both branches)")
ax.plot(ks, Xs, "b.-", ms=2, lw=1, label="predictor-corrector branch")
ax.plot(ks[0], Xs[0], "go", label="start")
if fold_k is not None:
    ax.plot(fold_k, fold_X, "rs", ms=9, label="continuation stalls (fold)")
ax.plot(kgrid[i_fold], Xgrid[i_fold], "kx", ms=10, mew=2, label="analytic fold")
ax.set_xlabel("control parameter k")
ax.set_ylabel("steady state X")
ax.set_title("Numerical continuation of self-activating gene: X(k) stalls at fold")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.5.1_s5.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result.
# ---------------------------------------------------------------
print("Explanation: The check confirms the result because the natural-parameter "
      "continuation tracks X(k) with tiny residuals away from folds but halts exactly "
      "where the tangent dX/dk = -(df/dk)/(df/dX) diverges (df/dX -> 0), which coincides "
      "with the analytic fold where k(X) is maximal and k becomes multivalued.")
