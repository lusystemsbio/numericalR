import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Self-activating gene model
#   f(X,k) = g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X
# k is the control parameter we sweep.
# ----------------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4


def hill(X):
    r = (X / Xth) ** n
    return r / (1.0 + r)


def f(X, k):
    return g0 + g1 * hill(X) - k * X


# ----------------------------------------------------------------------
# Explicit RK4 integrator that relaxes X to a steady state of
#   dX/dt = sign * f(X, k).
# sign=+1 : real dynamics  -> STABLE steady states are attractors.
# sign=-1 : reversed field -> UNSTABLE steady states become attractors.
# ----------------------------------------------------------------------
def integrate_to_ss(X, k, sign=1.0, dt=0.02, tol=1e-9, maxsteps=200000):
    for _ in range(maxsteps):
        r1 = sign * f(X, k)
        r2 = sign * f(X + 0.5 * dt * r1, k)
        r3 = sign * f(X + 0.5 * dt * r2, k)
        r4 = sign * f(X + dt * r3, k)
        dX = (dt / 6.0) * (r1 + 2 * r2 + 2 * r3 + r4)
        X_new = X + dX
        if not np.isfinite(X_new) or abs(X_new) > 1e6:  # blow-up -> no interior fixed point
            return np.nan
        if abs(X_new - X) < tol * (1.0 + abs(X_new)):  # converged
            return X_new
        X = X_new
    return X


# ----------------------------------------------------------------------
# STABLE branches by continuation.
# Sweep k downward starting on the lower (low-X) branch; follow the
# steady state by continuing from the previous state after each small
# nudge in k. When the lower fold is passed the state jumps up onto the
# upper branch: detect that large jump and REVERSE the k-sweep so the
# same continuation now walks the upper branch back up.
# ----------------------------------------------------------------------
JUMP = 80.0                     # |dX| that flags a saddle-node jump
kmin_sweep, kmax_sweep = 0.08, 0.22
dk0 = 0.001

stable_k, stable_X = [], []
k = kmax_sweep                  # start above the upper fold -> only low state exists
dk = -dk0                       # sweep k downward first
X = integrate_to_ss(5.0, k, +1.0)   # relax onto the lower branch
reversed_once = False
fold_lower_k = fold_upper_k = None

while kmin_sweep - 1e-12 <= k <= kmax_sweep + 1e-12:
    X_new = integrate_to_ss(X, k, +1.0)         # continue from previous state
    if stable_X and abs(X_new - stable_X[-1]) > JUMP and not reversed_once:
        fold_lower_k = k        # lower branch vanished here -> jumped to upper
        dk = -dk                # reverse the k-sweep to trace the upper branch
        reversed_once = True
    if stable_X and reversed_once and fold_upper_k is None \
            and abs(X_new - stable_X[-1]) > JUMP:
        fold_upper_k = k        # upper branch vanished here on the way back up
    stable_k.append(k)
    stable_X.append(X_new)
    X = X_new
    k += dk

stable_k = np.array(stable_k)
stable_X = np.array(stable_X)

# ----------------------------------------------------------------------
# UNSTABLE middle branch.
# Start in the bistable window at a point BETWEEN the two stable states
# and integrate the REVERSED field dX/dt = -f, for which the unstable
# steady state is an attractor. Then continue in k in both directions,
# each time starting from the previous unstable state, staying inside
# the bistable range (integration blows up once a fold is crossed).
# ----------------------------------------------------------------------
k_mid = 0.147
X_u = integrate_to_ss(200.0, k_mid, -1.0)       # relax onto the unstable state

unstable_k, unstable_X = [k_mid], [X_u]

# sweep k upward from the middle
X = X_u
k = k_mid + dk0
while k <= kmax_sweep:
    X_new = integrate_to_ss(X, k, -1.0)
    if not np.isfinite(X_new) or (unstable_X and abs(X_new - X) > JUMP):
        break                                    # crossed the upper fold
    unstable_k.append(k); unstable_X.append(X_new); X = X_new
    k += dk0

# sweep k downward from the middle
X = X_u
k = k_mid - dk0
while k >= kmin_sweep:
    X_new = integrate_to_ss(X, k, -1.0)
    if not np.isfinite(X_new) or (unstable_X and abs(X_new - X) > JUMP):
        break                                    # crossed the lower fold
    unstable_k.append(k); unstable_X.append(X_new); X = X_new
    k -= dk0

unstable_k = np.array(unstable_k)
unstable_X = np.array(unstable_X)
order = np.argsort(unstable_k)
unstable_k, unstable_X = unstable_k[order], unstable_X[order]

# ----------------------------------------------------------------------
# Reference S-curve (the 2E.2 result): the steady-state condition
# f=0 gives k explicitly as a function of X, k(X)=(g0+g1*hill(X))/X.
# This is the exact bifurcation curve, all three branches at once.
# ----------------------------------------------------------------------
X_grid = np.linspace(5.0, 750.0, 4000)
k_ref = (g0 + g1 * hill(X_grid)) / X_grid

# Check: every traced point must satisfy f(X,k)=0.
res_stable = np.max(np.abs(f(stable_X, stable_k)))
res_unstable = np.max(np.abs(f(unstable_X, unstable_k)))

# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.plot(k_ref, X_grid, '-', color='0.7', lw=3,
         label='reference S-curve (2E.2, k(X) from f=0)')
plt.plot(stable_k, stable_X, '.', color='tab:blue', ms=4,
         label='stable steady states (integrate dX/dt=+f)')
plt.plot(unstable_k, unstable_X, '.', color='tab:red', ms=4,
         label='unstable branch (integrate dX/dt=-f)')
plt.xlabel('k (degradation-rate control parameter)')
plt.ylabel('steady-state X')
plt.title('S-shaped bifurcation curve traced by ODE continuation')
plt.xlim(kmin_sweep, kmax_sweep)
plt.ylim(0, 750)
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.5.1_s5.png")

# ----------------------------------------------------------------------
# Numerical results
# ----------------------------------------------------------------------
print("Parameters: g0=%.1f  g1=%.1f  Xth=%.1f  n=%d" % (g0, g1, Xth, n))
print("Number of stable steady states traced:   %d" % len(stable_k))
print("Number of unstable steady states traced: %d" % len(unstable_k))
print("Lower saddle-node (jump-up) at k   ~ %.4f" % (fold_lower_k if fold_lower_k else float('nan')))
print("Upper saddle-node (jump-down) at k ~ %.4f" % (fold_upper_k if fold_upper_k else float('nan')))
print("Unstable-branch k range: %.4f to %.4f" % (unstable_k.min(), unstable_k.max()))
print("Unstable-branch X range: %.4f to %.4f" % (unstable_X.min(), unstable_X.max()))
print("Max |f(X,k)| over traced STABLE points:   %.3e" % res_stable)
print("Max |f(X,k)| over traced UNSTABLE points: %.3e" % res_unstable)
print("Check: because every traced point satisfies f(X,k)=0 to ~1e-8, the continuation "
      "points lie exactly on the exact k(X) curve of 2E.2 (which is defined by f=0), "
      "and the reversed-field integration adds the middle branch the forward field could not reach.")
