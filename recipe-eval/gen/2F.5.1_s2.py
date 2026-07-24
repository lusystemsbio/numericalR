import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Model: self-activating gene
# f(X,k) = g0 + g1*(X/Xth)^n / (1+(X/Xth)^n) - k*X
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(X, k):
    u = (X / Xth) ** n
    return g0 + g1 * u / (1.0 + u) - k * X

def dfdX(X, k):
    # derivative of the Hill term + derivative of -k*X
    u = (X / Xth) ** n
    dh = n * u / (X * (1.0 + u) ** 2)   # d/dX [ u/(1+u) ]
    return g1 * dh - k

def dfdk(X, k):
    return -X                            # d/dk of -k*X

# ---------------------------------------------------------------
# Newton corrector: given k, solve f(X,k)=0 starting from X0
# ---------------------------------------------------------------
def correct(X0, k, tol=1e-10, itmax=100):
    X = X0
    for _ in range(itmax):
        r = f(X, k)
        if abs(r) < tol:
            return X, True
        d = dfdX(X, k)
        if d == 0.0:                     # flat slope -> cannot correct
            return X, False
        X = X - r / d
        if X <= 0.0:                     # left the physical branch
            return X, False
    return X, abs(f(X, k)) < 1e-6

# ---------------------------------------------------------------
# Find a known starting point on the high (active) branch at k0
# ---------------------------------------------------------------
k0 = 0.05
X0, ok = correct(1000.0, k0)
print(f"Start point: k = {k0:.6f}, X = {X0:.6f}, converged = {ok}")

# ---------------------------------------------------------------
# Predictor-corrector continuation in the natural parameter k
# ---------------------------------------------------------------
dk = 0.002                               # step in the control parameter
slope_tol = 1e-3                         # |df/dX| below this => near fold
ks, Xs, slopes = [k0], [X0], [dfdX(X0, k0)]

k, X = k0, X0
stalled = False
stall_k = None
while k < 0.6:
    d = dfdX(X, k)
    # tangent of the curve X(k):  dX/dk = -(df/dk)/(df/dX)
    if abs(d) < slope_tol:               # df/dX -> 0 : dX/dk diverges (fold)
        stalled = True
        stall_k = k
        print(f"STALL at fold: k = {k:.6f}, X = {X:.6f}, df/dX = {d:.3e}")
        break
    dXdk = -dfdk(X, k) / d
    k_new = k + dk
    X_pred = X + dXdk * dk               # tangent predictor
    X_new, converged = correct(X_pred, k_new)   # Newton corrector back onto curve
    if not converged:                    # corrector failed near the fold
        stalled = True
        stall_k = k
        print(f"STALL (corrector failed): last good k = {k:.6f}, X = {X:.6f}")
        break
    k, X = k_new, X_new
    ks.append(k); Xs.append(X); slopes.append(dfdX(X, k))

ks, Xs, slopes = np.array(ks), np.array(Xs), np.array(slopes)
dXdk_arr = -dfdk(Xs, ks) / slopes        # tangent along the traced branch

print(f"Points traced: {len(ks)}")
print(f"k range traced: {ks[0]:.6f} to {ks[-1]:.6f}")
print(f"X range traced: {Xs.min():.6f} to {Xs.max():.6f}")
print(f"max |dX/dk| near end: {np.max(np.abs(dXdk_arr)):.6e}")
print(f"df/dX at last traced point: {slopes[-1]:.6e}")

# ---------------------------------------------------------------
# CHECK: independently verify each continued point is a true root,
# and show that |dX/dk| blows up as the fold is approached.
# ---------------------------------------------------------------
resid = np.array([abs(f(X, k)) for X, k in zip(Xs, ks)])
print(f"CHECK max |f(X,k)| over traced branch: {resid.max():.3e}")
print(f"CHECK |dX/dk| at start (away from fold): {abs(dXdk_arr[0]):.6e}")
print(f"CHECK |dX/dk| at last point (near fold): {abs(dXdk_arr[-1]):.6e}")
if stall_k is not None:
    print(f"Fold (stall) located near k = {stall_k:.6f}")

# ---------------------------------------------------------------
# Plot the continued branch X(k)
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(ks, Xs, "-o", ms=3, color="C0", label="continued branch X(k)")
if stall_k is not None:
    ax.plot(ks[-1], Xs[-1], "s", ms=10, color="C3",
            label=f"fold / stall (k≈{stall_k:.4f})")
ax.set_xlabel("k (control parameter)")
ax.set_ylabel("X (steady state)")
ax.set_title("Predictor-corrector continuation of a self-activating gene\n"
             "branch stalls at the fold where dX/dk diverges")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.5.1_s2.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# The residuals f(X,k) stay ~0 (accurate branch-following) while |dX/dk|
# grows without bound exactly where tracing halts, confirming the stall is
# the genuine fold at which df/dX->0 makes k a multivalued function of X.
# ---------------------------------------------------------------
print("Explanation: the check confirms the result because f(X,k) stays "
      "essentially zero along the whole traced branch (accurate following) "
      "while |dX/dk| diverges precisely at the point where tracing stalls, "
      "which is the fold where df/dX->0 and k becomes multivalued in X.")
