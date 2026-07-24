import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Model: self-activating gene
#   f(X,k) = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4

def f(X, k):
    """RHS of dX/dt = f(X,k): basal + Hill activation - linear degradation."""
    h = (X / Xth) ** n
    return g0 + g1 * h / (1.0 + h) - k * X

# ---------------------------------------------------------------
# Explicit relaxation-to-steady-state integrator (RK4, fixed step).
# sign=+1 -> integrate dX/dt=+f : STABLE steady states are attractors.
# sign=-1 -> integrate dX/dt=-f : time-reversed, so the UNSTABLE
#            (middle) steady state becomes an attractor and can be found.
# ---------------------------------------------------------------
def relax(X, k, sign, dt=0.05, nmax=400000, tol=1e-9, Xcap=3000.0):
    for _ in range(nmax):
        # one RK4 step of dX/dt = sign*f(X,k)
        k1 = sign * f(X, k)
        k2 = sign * f(X + 0.5 * dt * k1, k)
        k3 = sign * f(X + 0.5 * dt * k2, k)
        k4 = sign * f(X + dt * k3, k)
        dX = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0
        X = X + dt * dX
        if X < 0.0:              # concentration cannot go negative
            X = 0.0
        if X > Xcap:             # trajectory escaped (no attractor here)
            return X, False
        if abs(dX) < tol:        # converged to a steady state
            return X, True
    return X, True

# ---------------------------------------------------------------
# Sweep settings across the bistable range of k
# ---------------------------------------------------------------
kmin, kmax, dk = 0.10, 0.19, 0.0015
jump = 40.0   # a change in X larger than this = we fell off a fold

# ===============================================================
# STEP 1 -- LOWER STABLE BRANCH: start high k / low X, sweep k DOWN.
# Integrate +f to the stable low state, nudge k, continue from the
# previous state.  When X jumps up we have passed the LEFT fold.
# ===============================================================
lower = []
X, _ = relax(g0 / kmax, kmax, +1)      # settle onto the low state at k=kmax
k = kmax
k_fold_left, X_fold_left = None, None
while k >= kmin - 1e-12:
    Xnew, _ = relax(X, k, +1)          # relax to stable state, continuing from X
    if lower and abs(Xnew - lower[-1][1]) > jump:   # big jump -> left fold reached
        k_fold_left, X_fold_left = lower[-1]
        break
    lower.append((k, Xnew))
    X = Xnew
    k -= dk                            # nudge the control parameter

# ===============================================================
# STEP 2 -- REVERSE the sweep at the fold and switch to dX/dt=-f.
# Starting from the near-fold state (where lower & middle branches
# merge), integrate -f while sweeping k UP: the unstable middle
# branch is now an attractor, so we trace it until it merges with
# the upper branch at the RIGHT fold (trajectory then escapes).
# ===============================================================
middle = []
X = X_fold_left
k = k_fold_left
k_fold_right = None
while k <= kmax + 1e-12:
    Xnew, ok = relax(X, k, -1)         # -f makes the unstable state attracting
    if not ok:                         # escaped -> past the right fold
        k_fold_right = middle[-1][0] if middle else k
        break
    if middle and abs(Xnew - middle[-1][1]) > jump:
        k_fold_right = middle[-1][0]
        break
    middle.append((k, Xnew))
    X = Xnew
    k += dk                            # nudge k back the other way

# ===============================================================
# STEP 3 -- UPPER STABLE BRANCH: start low k / high X, sweep k UP.
# Integrate +f to the stable high state, continue from previous X,
# until X jumps down (RIGHT fold), for completeness of the S-curve.
# ===============================================================
upper = []
X, _ = relax((g0 + g1) / kmin, kmin, +1)   # settle onto the high state at k=kmin
k = kmin
while k <= kmax + 1e-12:
    Xnew, _ = relax(X, k, +1)
    if upper and abs(Xnew - upper[-1][1]) > jump:   # jump down -> right fold
        if k_fold_right is None:
            k_fold_right = upper[-1][0]
        break
    upper.append((k, Xnew))
    X = Xnew
    k += dk

lower = np.array(lower)
middle = np.array(middle)
upper = np.array(upper)

# ===============================================================
# CHECK vs 2E.2: the exact steady-state locus is k = f-balance,
# i.e. k(X) = (g0 + g1*H(X))/X, obtained by setting f(X,k)=0.
# The traced points must satisfy f(X,k)=0 -> lie on this S-curve.
# ===============================================================
Xg = np.linspace(1.0, 800.0, 4000)
hg = (Xg / Xth) ** n
kg = (g0 + g1 * hg / (1.0 + hg)) / Xg      # exact S-curve of 2E.2

all_pts = np.vstack([lower, middle, upper])
residuals = np.array([abs(f(Xp, kp)) for kp, Xp in all_pts])
max_res = residuals.max()

# ---------------------------------------------------------------
# Print numerical results
# ---------------------------------------------------------------
print(f"Left (saddle-node) fold:  k = {k_fold_left:.5f},  X = {X_fold_left:.3f}")
print(f"Right (saddle-node) fold: k = {k_fold_right:.5f}")
print(f"Bistable range of k: [{k_fold_left:.5f}, {k_fold_right:.5f}]")
print(f"Lower stable branch points:   {len(lower)}")
print(f"Middle unstable branch points:{len(middle)}")
print(f"Upper stable branch points:   {len(upper)}")
print(f"Lower branch X range: [{lower[:,1].min():.3f}, {lower[:,1].max():.3f}]")
print(f"Middle branch X range:[{middle[:,1].min():.3f}, {middle[:,1].max():.3f}]")
print(f"Upper branch X range: [{upper[:,1].min():.3f}, {upper[:,1].max():.3f}]")
print(f"Max |f(X,k)| over ALL traced points (should be ~0): {max_res:.3e}")

# max distance of traced points from the exact 2E.2 locus (in k, at fixed X)
def k_exact(X):
    h = (X / Xth) ** n
    return (g0 + g1 * h / (1.0 + h)) / X
max_k_dev = max(abs(kp - k_exact(Xp)) for kp, Xp in all_pts)
print(f"Max |k_traced - k_exact(X)| vs 2E.2 curve: {max_k_dev:.3e}")

# ---------------------------------------------------------------
# Plot
# ---------------------------------------------------------------
out = "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.5.1_s1.png"
os.makedirs(os.path.dirname(out), exist_ok=True)

plt.figure(figsize=(7, 5))
plt.plot(kg, Xg, color="0.8", lw=6, label="exact S-curve (2E.2)")
plt.plot(lower[:, 0], lower[:, 1], "b-", lw=2, label="lower stable branch")
plt.plot(upper[:, 0], upper[:, 1], "b-", lw=2)
plt.plot(middle[:, 0], middle[:, 1], "r--", lw=2, label="unstable middle branch (-f)")
plt.scatter([k_fold_left, k_fold_right],
            [X_fold_left, np.interp(k_fold_right, upper[:, 0], upper[:, 1])],
            color="k", zorder=5, label="folds")
plt.xlim(kmin, kmax)
plt.ylim(0, 700)
plt.xlabel("k (control parameter)")
plt.ylabel("steady-state X")
plt.title("S-shaped bifurcation traced by ODE integration")
plt.legend(loc="upper right")
plt.tight_layout()
plt.savefig(out)

# One-sentence explanation of why the check confirms the result:
print("Explanation: every traced point has |f(X,k)|~0, so each lies exactly on "
      "the steady-state locus k=(g0+g1*H(X))/X plotted in 2E.2, and the -f-traced "
      "points fill the previously-missing unstable middle segment, reproducing the "
      "full S-curve including its unstable branch.")
