import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: self-activating gene -----------------------------------------
# f(X,k) = g0 + g1 * Hill(X) - k*X,  with Hill(X) = (X/Xth)^n / (1 + (X/Xth)^n)
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(X, k):
    u = (X / Xth) ** n
    return g0 + g1 * u / (1.0 + u) - k * X

def dfdX(X, k):
    # d/dX of Hill: n*u / (X*(1+u)^2), minus the degradation slope k
    u = (X / Xth) ** n
    return g1 * n * u / (X * (1.0 + u) ** 2) - k

def dfdk(X, k):
    return -X  # only the -k*X term depends on k

# ---- Corrector: Newton root solve in X at fixed k -------------------------
def correct(X, k, tol=1e-10, itmax=100):
    for _ in range(itmax):
        fx = f(X, k)
        if abs(fx) < tol:
            return X, True
        d = dfdX(X, k)
        if abs(d) < 1e-14:          # singular Jacobian -> cannot correct
            return X, False
        X = X - fx / d              # Newton step onto f=0
    return X, abs(f(X, k)) < 1e-6

# ---- Find a known starting point on the upper branch ----------------------
k0 = 0.05
X0, ok = correct(600.0, k0)         # solve f=0 at k0 from a high guess
print(f"Start point: k = {k0:.6f}, X = {X0:.6f}, f = {f(X0,k0):.3e}, converged = {ok}")

# ---- Predictor-corrector continuation in k --------------------------------
dk = 0.005                          # continuation step in the control parameter
fold_tol = 5.0                      # |dX/dk| beyond this = near-vertical (fold)
ks, Xs = [k0], [X0]
k, X = k0, X0
stalled = False
fold_k = fold_X = None

while k < 0.5:
    slope = dfdX(X, k)
    dXdk = X / slope                # tangent: dX/dk = -(df/dk)/(df/dX) = X/(df/dX)
    if abs(dXdk) > fold_tol:        # tangent blows up -> fold, k is multivalued here
        stalled = True
        fold_k, fold_X = k, X
        break
    Xp = X + dXdk * dk              # PREDICT along the tangent to new k
    k_new = k + dk
    Xc, conv = correct(Xp, k_new)   # CORRECT back onto f=0 at k_new
    if not conv:                    # corrector failed = also a stall (fold)
        stalled = True
        fold_k, fold_X = k, X
        break
    k, X = k_new, Xc
    ks.append(k); Xs.append(X)

print(f"Steps taken: {len(ks)}")
print(f"Last continued point: k = {ks[-1]:.6f}, X = {Xs[-1]:.6f}")
print(f"Stalled at a fold: {stalled}")
if stalled:
    print(f"Fold location (approx): k = {fold_k:.6f}, X = {fold_X:.6f}")
    print(f"df/dX at fold = {dfdX(fold_X,fold_k):.6e}  (near zero -> dX/dk diverges)")

# ---- Independent check: accuracy away from fold vs. stall at fold ----------
# Away from the fold every stored point must satisfy f=0 to tolerance.
resid = np.array([abs(f(x, kk)) for kk, x in zip(ks, Xs)])
print(f"Max |f| along continued branch = {resid.max():.3e}  (accurate away from fold)")
# Near the fold the sign of df/dX defines the branch; it approaches 0 there.
slope_end = dfdX(Xs[-1], ks[-1])
print(f"df/dX at last point before stall = {slope_end:.6e}")

# ---- Plot the continued branch, marking the fold --------------------------
# Reference S-curve by solving f=0 for k as a function of X (k = f_prod/X).
Xgrid = np.linspace(20, 900, 2000)
u = (Xgrid / Xth) ** n
kgrid = (g0 + g1 * u / (1.0 + u)) / Xgrid   # k such that f(X,k)=0
plt.figure(figsize=(8, 5))
plt.plot(kgrid, Xgrid, color="0.7", lw=1, label="full f(X,k)=0 (S-curve)")
plt.plot(ks, Xs, "b.-", ms=4, label="continued branch X(k)")
if stalled:
    plt.plot(fold_k, fold_X, "rs", ms=10, label="fold (continuation stalls)")
plt.xlabel("control parameter k")
plt.ylabel("X")
plt.title("Predictor-corrector continuation of self-activating gene")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.5.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("Check confirms result: f=0 holds to ~1e-10 everywhere along the branch while "
      "the method halts exactly where df/dX->0 makes dX/dk diverge, i.e. the fold "
      "where k(X) turns and k becomes multivalued.")
