import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Self-activating gene: basal + excitatory Hill - linear degradation
# f(X,k) = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X ;  k is the control parameter
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4

def f(X, k):
    r = (X / Xth) ** n
    return g0 + g1 * r / (1.0 + r) - k * X

def k_of_X(X):
    # Exact steady-state relation from f(X,k)=0 :  k = (g0 + g1 H(X)) / X
    r = (X / Xth) ** n
    return (g0 + g1 * r / (1.0 + r)) / X

# --- explicit RK4 relaxation of dX/dt = sign*f(X,k) to a steady state ---
# sign=+1 : flow of the real ODE  -> stable steady states are attractors
# sign=-1 : flow of the time-reversed ODE -> the UNSTABLE steady state
#           (df/dX>0) becomes an attractor, so we can "sit on" it.
def relax(X0, k, sign=+1, dt=0.1, tol=1e-6, max_steps=500000):
    X = float(X0)
    for _ in range(max_steps):
        F1 = sign * f(X, k)
        if abs(F1) < tol:                    # reached a steady state
            break
        F2 = sign * f(X + 0.5 * dt * F1, k)
        F3 = sign * f(X + 0.5 * dt * F2, k)
        F4 = sign * f(X + dt * F3, k)
        X += (dt / 6.0) * (F1 + 2 * F2 + 2 * F3 + F4)
    return X

# ---------------------------------------------------------------
# Locate the two saddle-node folds (bounds of the bistable range)
# as the vertical-tangent points of the exact S-curve k(X).
# ---------------------------------------------------------------
Xgrid = np.linspace(1.0, 3000.0, 300000)
kgrid = k_of_X(Xgrid)
slope_sign = np.sign(np.diff(kgrid))
fold_idx = np.where(np.diff(slope_sign) != 0)[0] + 1
folds = np.sort(kgrid[fold_idx])
k_a, k_b = folds[0], folds[-1]          # lower fold, upper fold
print(f"Lower saddle-node fold  k_a = {k_a:.6f}")
print(f"Upper saddle-node fold  k_b = {k_b:.6f}")

# Sweep a little outside the folds so both monostable tails show up.
pad = 0.15 * (k_b - k_a)
k_min, k_max = k_a - pad, k_b + pad
print(f"Sweep range  k_min = {k_min:.6f}  k_max = {k_max:.6f}")

# ---------------------------------------------------------------
# 1) STABLE BRANCHES by ODE integration with reversal at a jump.
#    Start high on the upper branch, march k up until the upper
#    branch vanishes (large jump down), reverse the sweep, march k
#    back down along the lower branch until it too jumps.
# ---------------------------------------------------------------
dk = (k_max - k_min) / 200.0
jump = 80.0                                   # "large jump" detector

k = k_min
X = relax(1500.0, k, sign=+1)                 # settle onto the upper branch
step = +dk
reversals = 0
ks_stab, Xs_stab = [], []
while True:
    ks_stab.append(k); Xs_stab.append(X)      # record current steady state
    kn = k + step                             # nudge the control parameter
    if kn < k_min or kn > k_max:              # hit a sweep boundary
        break
    Xn = relax(X, kn, sign=+1)                # continue from the previous state
    if abs(Xn - X) > jump:                    # a branch disappeared -> jump
        reversals += 1
        step = -step                          # reverse the k-sweep direction
        k, X = kn, Xn                         # keep the post-jump point
        if reversals >= 2:                    # both folds seen -> done
            ks_stab.append(k); Xs_stab.append(X)
            break
        continue
    k, X = kn, Xn
ks_stab = np.array(ks_stab); Xs_stab = np.array(Xs_stab)

# ---------------------------------------------------------------
# 2) UNSTABLE MIDDLE BRANCH: integrate the time-reversed ODE
#    dX/dt = -f, for which the middle (df/dX>0) state is an attractor.
#    Seed inside the bistable window, then continue in k both ways.
# ---------------------------------------------------------------
eps = 0.02 * (k_b - k_a)
k_mid = 0.5 * (k_a + k_b)
X_mid = relax(Xth, k_mid, sign=-1)            # sit on the unstable state via -f
ks_uns, Xs_uns = [k_mid], [X_mid]

Xc = X_mid                                    # continue upward in k
kc = k_mid
while kc + dk < k_b - eps:
    kc += dk
    Xc = relax(Xc, kc, sign=-1)               # follow unstable state as attractor of -f
    ks_uns.append(kc); Xs_uns.append(Xc)

Xc = X_mid                                    # continue downward in k
kc = k_mid
while kc - dk > k_a + eps:
    kc -= dk
    Xc = relax(Xc, kc, sign=-1)
    ks_uns.append(kc); Xs_uns.append(Xc)

order = np.argsort(ks_uns)
ks_uns = np.array(ks_uns)[order]; Xs_uns = np.array(Xs_uns)[order]

print(f"Stable-branch points traced   = {len(ks_stab)}")
print(f"Unstable-branch points traced = {len(ks_uns)}")

# ---------------------------------------------------------------
# 3) CHECK: every traced (k,X) must satisfy f(X,k)=0, i.e. lie on the
#    exact S-curve of 2E.2  ( k = (g0+g1 H(X))/X ).
# ---------------------------------------------------------------
res_stab = np.abs(f(Xs_stab, ks_stab))
res_uns = np.abs(f(Xs_uns, ks_uns))
print(f"Max |f(X,k)| on stable branches   = {res_stab.max():.3e}")
print(f"Max |f(X,k)| on unstable branch   = {res_uns.max():.3e}")
print(f"Max |f(X,k)| over ALL traced pts  = {max(res_stab.max(), res_uns.max()):.3e}")

# ---------------------------------------------------------------
# Plot: traced points over the exact analytic S-curve.
# ---------------------------------------------------------------
Xline = np.linspace(3.0, 2600.0, 4000)
kline = k_of_X(Xline)
inwin = (kline >= k_min) & (kline <= k_max)

plt.figure(figsize=(7, 5))
plt.plot(kline[inwin], Xline[inwin], '-', color='0.6', lw=1.0,
         label='exact S-curve (2E.2)')
plt.plot(ks_stab, Xs_stab, 'o', ms=3, color='tab:blue', label='traced stable (dX/dt=+f)')
plt.plot(ks_uns, Xs_uns, 's', ms=3, color='tab:red', label='traced unstable (dX/dt=-f)')
plt.axvline(k_a, ls=':', color='0.4'); plt.axvline(k_b, ls=':', color='0.4')
plt.xlabel('k (control parameter)'); plt.ylabel('steady-state X')
plt.title('S-shaped bifurcation traced by ODE integration')
plt.legend(loc='upper right', fontsize=8)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.5.1_s4.png")

# One-sentence explanation of the check:
print("Check rationale: since |f(X,k)|~0 at every traced point, each point is an "
      "exact steady state lying on the algebraic curve k=(g0+g1 H(X))/X of 2E.2, so "
      "recovering all three branches (upper, lower, and the middle branch caught via -f) "
      "confirms the continuation reproduced the full S-curve including its unstable arm.")
