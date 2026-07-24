import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model: self-activating gene -----
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4

def hill(X):
    u = (X / Xth)**n
    return u / (1.0 + u)

def f(X, k):
    # basal + excitatory Hill - linear degradation
    return g0 + g1 * hill(X) - k * X

def dfdX(X, k):
    # d/dX of Hill: u=(X/Xth)^n, dH/dX = (n*u/X)/(1+u)^2
    u = (X / Xth)**n
    dH = (n * u / X) / (1.0 + u)**2
    return g1 * dH - k

def dfdk(X, k):
    return -X

# ----- Corrector: Newton root solve for f(X,k)=0 at fixed k -----
def correct(X0, k, tol=1e-10, itmax=100):
    X = X0
    for _ in range(itmax):
        fx = f(X, k)
        dfx = dfdX(X, k)
        if abs(dfx) < 1e-14:          # Jacobian ~ 0 => at/near fold, cannot correct
            return X, False
        step = fx / dfx
        X -= step
        if X <= 0:                    # keep on the physical (positive) branch
            return X, False
        if abs(step) < tol:
            return X, True
    return X, False

# ----- Find a starting point on the upper branch at a small k -----
k0 = 0.15
# scan for the largest positive root at k0
Xs = np.linspace(1.0, 1200.0, 20000)
fv = f(Xs, k0)
sign = np.sign(fv)
roots = []
for i in range(len(Xs) - 1):
    if sign[i] * sign[i + 1] < 0:
        Xr, ok = correct(0.5 * (Xs[i] + Xs[i + 1]), k0)
        if ok:
            roots.append(Xr)
X0 = max(roots)   # upper branch
print(f"Start point: k = {k0:.6f}, X = {X0:.6f}")
print(f"Residual f(X0,k0) = {f(X0, k0):.3e}")

# ----- Predictor-corrector continuation in k (increasing k toward the fold) -----
dk = 0.001
k_curr, X_curr = k0, X0
ks, Xk, dXdk_list = [k_curr], [X_curr], []
stall_k = None

while k_curr < 1.0:
    dfx = dfdX(X_curr, k_curr)
    # tangent slope dX/dk = -(df/dk)/(df/dX)
    if abs(dfx) < 1e-6:              # slope diverges: fold reached, k is multivalued -> stall
        stall_k = k_curr
        break
    slope = -dfdk(X_curr, k_curr) / dfx
    dXdk_list.append(slope)
    # Predict along the tangent
    k_next = k_curr + dk
    X_pred = X_curr + slope * dk
    # Correct back onto f=0 at the new k
    X_new, ok = correct(X_pred, k_next)
    if (not ok) or X_new <= 0:       # corrector cannot land on the curve near the fold
        stall_k = k_curr
        break
    k_curr, X_curr = k_next, X_new
    ks.append(k_curr)
    Xk.append(X_curr)

print(f"Points computed on branch: {len(ks)}")
print(f"Stalled at k (fold) = {stall_k:.6f}")
print(f"Last X before stall  = {Xk[-1]:.6f}")
print(f"Max |dX/dk| observed = {max(abs(np.array(dXdk_list))):.6e}")

# ----- Separate check: accuracy away from fold, divergence of dX/dk at fold -----
# Pick a point well away from the fold and compare tangent slope to a finite difference.
i_mid = len(ks) // 3
k_chk, X_chk = ks[i_mid], Xk[i_mid]
slope_tan = -dfdk(X_chk, k_chk) / dfdX(X_chk, k_chk)
h = 1e-4
Xp, _ = correct(X_chk, k_chk + h)
Xm, _ = correct(X_chk, k_chk - h)
slope_fd = (Xp - Xm) / (2 * h)
print(f"Check away from fold: k = {k_chk:.6f}")
print(f"  tangent dX/dk        = {slope_tan:.6f}")
print(f"  finite-diff dX/dk    = {slope_fd:.6f}")
print(f"  agreement error      = {abs(slope_tan - slope_fd):.3e}")
print(f"df/dX at last point (near fold) = {dfdX(Xk[-1], ks[-1]):.6e}")
print(f"|dX/dk| at last point (near fold) = {abs(slope_tan) if False else abs(-dfdk(Xk[-1], ks[-1])/dfdX(Xk[-1], ks[-1])):.6e}")

# ----- Plot -----
plt.figure(figsize=(8, 5))
plt.plot(ks, Xk, '-', color='C0', lw=2, label='continued branch X(k)')
plt.plot(ks[0], Xk[0], 'o', color='green', label='start')
plt.plot(ks[-1], Xk[-1], 's', color='red', label='stall at fold')
if stall_k is not None:
    plt.axvline(stall_k, ls='--', color='red', alpha=0.5,
                label=f'fold k ~ {stall_k:.3f}')
plt.xlabel('k (control parameter)')
plt.ylabel('X (steady state)')
plt.title('Predictor-corrector continuation of self-activating gene fixed points')
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2F.5.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The tangent slope dX/dk matches the finite-difference slope to high "
      "accuracy on the smooth part of the branch but blows up as df/dX -> 0 at the fold, "
      "confirming the continuation is correct away from folds and that stalling occurs "
      "exactly where k becomes multivalued and dX/dk diverges.")
