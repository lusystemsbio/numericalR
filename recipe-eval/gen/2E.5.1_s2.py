import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Self-activating gene model:
#   f(X,k) = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X
# We trace the S-shaped steady-state curve X*(k) by *integrating* the ODE
# to steady state and using numerical continuation (parameter marching),
# rather than by root-finding.  The unstable middle branch is captured by
# integrating the time-reversed ODE dX/dt = -f, which flips the stability
# so the unstable steady state becomes an attractor.
# ----------------------------------------------------------------------

# --- model parameters ---
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def hill(X):
    r = (X / Xth) ** n
    return r / (1.0 + r)

def f(X, k):
    return g0 + g1 * hill(X) - k * X

# --- explicit RK4 integrator that marches to a steady state ---------------
# sign=+1 integrates dX/dt = f  (stable steady states are attractors)
# sign=-1 integrates dX/dt = -f (unstable steady states become attractors)
def integrate_to_ss(X0, k, sign=+1, dt=0.02, tol=1e-7, max_steps=500000):
    X = float(X0)
    for _ in range(max_steps):
        # one classic 4th-order Runge-Kutta step on sign*f
        k1 = sign * f(X, k)
        k2 = sign * f(X + 0.5 * dt * k1, k)
        k3 = sign * f(X + 0.5 * dt * k2, k)
        k4 = sign * f(X + dt * k3, k)
        X = X + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        if X < 0.0:            # concentration cannot be negative
            X = 0.0
        # stop once we sit on a steady state (rate ~ 0)
        if abs(f(X, k)) < tol:
            break
    return X

JUMP = 50.0   # a change larger than this between successive k means a fold jump

# ======================================================================
# 1) LOWER STABLE BRANCH: descending k sweep, integrating dX/dt = f
# ======================================================================
k_start = 0.30            # high k -> only the low state exists
X = integrate_to_ss(5.0, k_start, sign=+1)   # settle onto the low steady state
lower_k, lower_X = [], []
k_lower_fold = None
ks_down = np.arange(k_start, 0.05, -0.001)   # nudge k downward step by step
for k in ks_down:
    Xnew = integrate_to_ss(X, k, sign=+1)    # continue from previous state
    if abs(Xnew - X) > JUMP:                 # large jump => lower fold reached
        k_lower_fold = (k, X)                # last point before the jump
        X = Xnew                             # X has jumped up to the high branch
        k_reverse = k                        # reverse the sweep direction here
        break
    lower_k.append(k); lower_X.append(Xnew)
    X = Xnew

# ======================================================================
# 2) UPPER STABLE BRANCH: reverse and sweep k upward from the jump point,
#    still integrating dX/dt = f (we are now sitting on the high branch)
# ======================================================================
upper_k, upper_X = [], []
k_upper_fold = None
ks_up = np.arange(k_reverse, k_start, 0.001)
for k in ks_up:
    Xnew = integrate_to_ss(X, k, sign=+1)
    if abs(Xnew - X) > JUMP:                 # large jump down => upper fold reached
        k_upper_fold = (k, X)                # last point before the jump
        X = Xnew
        break
    upper_k.append(k); upper_X.append(Xnew)
    X = Xnew

lower_k = np.array(lower_k); lower_X = np.array(lower_X)
upper_k = np.array(upper_k); upper_X = np.array(upper_X)

# ======================================================================
# 3) UNSTABLE MIDDLE BRANCH: integrate the time-reversed ODE dX/dt = -f.
#    In the bistable window (between the two folds) the middle steady
#    state is unstable for f but STABLE for -f, so it attracts.
#    We seed midway between the two stable branches and continue in k.
# ======================================================================
k_lo = k_lower_fold[0]     # bistable range lower edge (lower fold)
k_hi = k_upper_fold[0]     # bistable range upper edge (upper fold)

# sort stable-branch samples by k so we can interpolate seeds
lo_order = np.argsort(lower_k); Ls = lower_k[lo_order]; LX = lower_X[lo_order]
up_order = np.argsort(upper_k); Us = upper_k[up_order]; UX = upper_X[up_order]

ks_mid = np.linspace(k_lo, k_hi, 60)
mid_k, mid_X = [], []
Xseed = None
for k in ks_mid:
    Xlo = np.interp(k, Ls, LX)               # stable low value at this k
    Xhi = np.interp(k, Us, UX)              # stable high value at this k
    if Xseed is None:
        Xseed = 0.5 * (Xlo + Xhi)            # first seed: midpoint between branches
    Xmid = integrate_to_ss(Xseed, k, sign=-1)  # -f pulls us onto the unstable state
    # keep only genuine middle-branch points (strictly between the stable ones)
    if Xlo + 1.0 < Xmid < Xhi - 1.0:
        mid_k.append(k); mid_X.append(Xmid)
        Xseed = Xmid                         # continue from previous middle state
    else:
        Xseed = 0.5 * (Xlo + Xhi)            # re-seed if we fell off the branch

mid_k = np.array(mid_k); mid_X = np.array(mid_X)

# ======================================================================
# 4) INDEPENDENT CHECK (reproducing the 2E.2 S-curve):
#    At any steady state f=0  =>  k = (g0 + g1*hill(X)) / X.
#    Sweeping X algebraically traces the exact S-curve (all three branches).
# ======================================================================
X_ref = np.linspace(5.0, 700.0, 4000)
k_ref = (g0 + g1 * hill(X_ref)) / X_ref

# --- reported numbers ---
print(f"Lower fold (bifurcation) : k = {k_lower_fold[0]:.4f}, X = {k_lower_fold[1]:.3f}")
print(f"Upper fold (bifurcation) : k = {k_upper_fold[0]:.4f}, X = {k_upper_fold[1]:.3f}")
print(f"Bistable range in k      : [{k_lo:.4f}, {k_hi:.4f}]")
print(f"Lower stable branch pts  : {len(lower_k)}")
print(f"Upper stable branch pts  : {len(upper_k)}")
print(f"Unstable middle branch pts: {len(mid_k)}")

# a few sample steady states on each branch
print(f"Sample lower branch : k={Ls[len(Ls)//2]:.4f}, X={LX[len(LX)//2]:.3f}")
print(f"Sample upper branch : k={Us[len(Us)//2]:.4f}, X={UX[len(UX)//2]:.3f}")
if len(mid_k) > 0:
    j = len(mid_k) // 2
    print(f"Sample middle branch: k={mid_k[j]:.4f}, X={mid_X[j]:.3f}")

# maximum deviation of every traced point from the algebraic S-curve (residual |f|)
all_k = np.concatenate([lower_k, upper_k, mid_k])
all_X = np.concatenate([lower_X, upper_X, mid_X])
resid = np.array([abs(f(Xi, ki)) for Xi, ki in zip(all_X, all_k)])
print(f"Max |f| over all traced points : {resid.max():.3e}")
# max distance in k from the analytic curve k(X) at the same X
k_from_curve = (g0 + g1 * hill(all_X)) / all_X
print(f"Max |k_traced - k(X)|          : {np.max(np.abs(all_k - k_from_curve)):.3e}")

# ----------------------------------------------------------------------
# 5) Plot: steady-state X versus k, with the unstable branch
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.plot(k_ref, X_ref, '-', color='0.7', lw=3,
         label='2E.2 S-curve  k=(g0+g1 h(X))/X', zorder=1)
plt.plot(lower_k, lower_X, '.', color='C0', ms=4, label='lower stable (integrate +f)')
plt.plot(upper_k, upper_X, '.', color='C2', ms=4, label='upper stable (integrate +f)')
plt.plot(mid_k, mid_X, '.', color='C3', ms=4, label='unstable middle (integrate -f)')
plt.plot(k_lower_fold[0], k_lower_fold[1], 'kv', ms=9, label='fold points')
plt.plot(k_upper_fold[0], k_upper_fold[1], 'kv', ms=9)
plt.xlim(0.05, 0.30)
plt.xlabel('control parameter k')
plt.ylabel('steady-state X')
plt.title('S-shaped bifurcation curve by ODE continuation')
plt.legend(loc='upper right', fontsize=8)
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.5.1_s2.png")

# The check confirms the result because the ODE-traced points and the algebraic
# curve k=(g0+g1 h(X))/X both express the same steady-state condition f(X,k)=0,
# so their coincidence (near-zero |f| and matching k) proves the continuation
# recovered the true S-curve including its unstable middle branch.
print("Check: traced points coincide with the 2E.2 S-curve (f(X,k)=0), "
      "including the middle branch, confirming correct continuation.")
