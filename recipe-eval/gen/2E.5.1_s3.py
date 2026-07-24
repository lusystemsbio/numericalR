import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Self-activating gene model:
#   f(X,k) = basal + Hill-activation - linear degradation
# k (degradation rate) is the control parameter we sweep.
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4.0

def f(X, k):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X

def relax(X, k, sign, dt=0.05, tol=1e-8, nmax=2_000_000):
    """Integrate dX/dt = sign*f(X,k) with RK4 until velocity ~ 0.
    sign=+1 lands on a STABLE steady state of f;
    sign=-1 flips the flow so an UNSTABLE steady state of f becomes the attractor."""
    for _ in range(nmax):
        s1 = sign * f(X, k)                       # RK4 slopes
        s2 = sign * f(X + 0.5 * dt * s1, k)
        s3 = sign * f(X + 0.5 * dt * s2, k)
        s4 = sign * f(X + dt * s3, k)
        dX = (dt / 6.0) * (s1 + 2 * s2 + 2 * s3 + s4)
        X += dX
        if X < 0.0:                               # concentration stays physical
            X = 0.0
        if abs(dX / dt) < tol:                    # dX/dt ~ 0  -> steady state
            break
    return X

# ---------------------------------------------------------------
# PHASE 1: follow the UPPER stable branch while nudging k upward.
# Start from a high state, relax to steady, step k, continue from
# the previous state. A large downward jump = the upper fold.
# ---------------------------------------------------------------
upper = []
X = 500.0
reverse_k = None
prev = None
for k in np.arange(0.08, 0.22, 0.002):
    X = relax(X, k, +1.0)                         # integrate dX/dt = f to stable SS
    if prev is not None and X < prev - 40.0:      # big jump -> we fell off the fold
        reverse_k = k
        break
    upper.append((k, X))
    prev = X

# ---------------------------------------------------------------
# PHASE 2: reverse the sweep. We are now on the LOWER stable branch
# (the state jumped down). Step k downward, continuing each time,
# until a large upward jump = the lower fold.
# ---------------------------------------------------------------
lower = []
X = relax(X, reverse_k, +1.0)                     # settle onto the lower branch
prev = None
lower_fold_k = None
for k in np.arange(reverse_k, 0.08, -0.002):
    X = relax(X, k, +1.0)
    if prev is not None and X > prev + 40.0:      # big jump up -> lower fold reached
        lower_fold_k = k
        break
    lower.append((k, X))
    prev = X

# ---------------------------------------------------------------
# PHASE 3: capture the UNSTABLE middle branch. For each k in the
# bistable range, find both stable states, then integrate the
# REVERSED flow dX/dt = -f starting between them: the unstable
# steady state is now an attractor, so we can trace it.
# ---------------------------------------------------------------
middle = []
for k in np.arange(0.126, 0.170, 0.002):
    xl = relax(1.0,    k, +1.0)                   # lower stable state
    xu = relax(1000.0, k, +1.0)                   # upper stable state
    if abs(xu - xl) < 1.0:                        # not bistable at this k -> skip
        continue
    xm = relax(0.5 * (xl + xu), k, -1.0)          # -f flow -> unstable middle
    middle.append((k, xm))

# ---------------------------------------------------------------
# CHECK against the analytic S-curve of 2E.2:  f(X,k)=0  <=>
#   k(X) = (g0 + g1*Hill(X)) / X   (single-valued in X, S-shaped in k).
# ---------------------------------------------------------------
Xgrid = np.linspace(1.0, 800.0, 4000)
hill = (Xgrid / Xth) ** n / (1.0 + (Xgrid / Xth) ** n)
k_curve = (g0 + g1 * hill) / Xgrid

# max residual |f| over all traced points and max mismatch vs analytic k(X)
all_pts = upper + lower + middle
max_res = max(abs(f(X, k)) for (k, X) in all_pts)
max_kdev = max(abs(k - (g0 + g1 * (X / Xth) ** n / (1.0 + (X / Xth) ** n)) / X)
               for (k, X) in all_pts)

# ---------------------------------------------------------------
# Report numbers
# ---------------------------------------------------------------
print(f"Upper fold (jump-down) k          : {reverse_k:.4f}")
print(f"Lower fold (jump-up)   k          : {lower_fold_k:.4f}")
print(f"Bistable k-range approx           : [{lower_fold_k:.4f}, {reverse_k:.4f}]")
print(f"Upper-branch points traced        : {len(upper)}")
print(f"Lower-branch points traced        : {len(lower)}")
print(f"Middle (unstable) points traced   : {len(middle)}")
print(f"Sample upper SS  (k={upper[len(upper)//2][0]:.3f}) X = {upper[len(upper)//2][1]:.4f}")
print(f"Sample lower SS  (k={lower[len(lower)//2][0]:.3f}) X = {lower[len(lower)//2][1]:.4f}")
print(f"Sample middle SS (k={middle[len(middle)//2][0]:.3f}) X = {middle[len(middle)//2][1]:.4f}")
print(f"Max |f(X,k)| over traced points   : {max_res:.3e}")
print(f"Max |k_swept - k_analytic(X)|     : {max_kdev:.3e}")

# ---------------------------------------------------------------
# Plot
# ---------------------------------------------------------------
plt.figure(figsize=(7, 5))
plt.plot(k_curve, Xgrid, '-', color='0.6', lw=1.0, label='analytic S-curve (2E.2)')
ku, Xu = zip(*upper);   plt.plot(ku, Xu, 'o', ms=4, color='tab:blue',  label='upper stable (f)')
kl, Xl = zip(*lower);   plt.plot(kl, Xl, 's', ms=4, color='tab:green', label='lower stable (f)')
km, Xm = zip(*sorted(middle)); plt.plot(km, Xm, '^', ms=5, color='tab:red', label='unstable middle (-f)')
plt.xlabel('degradation rate  k')
plt.ylabel('steady-state  X')
plt.title('S-shaped bifurcation curve traced by ODE continuation')
plt.xlim(0.08, 0.22)
plt.ylim(0, 750)
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.5.1_s3.png")

# The check confirms the result because the analytic curve is exactly the locus f(X,k)=0
# (all steady states), so the near-zero residual |f| shows every traced point—including the
# reversed-flow middle branch—lands on the same S-curve, now completed with its unstable arm.
print("Check: traced points satisfy f(X,k)=0 (max residual above ~machine-small),")
print("so they lie on the 2E.2 S-curve and the -f pass fills in the unstable middle branch.")
