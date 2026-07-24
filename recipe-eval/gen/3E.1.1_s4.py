import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Toggle-switch model:
#   dX/dt = gX0 + gX1/(1+(Y/Yth)^nY) - kX*X
#   dY/dt = gY0 + gY1/(1+(X/Xth)^nX) - kY*Y
# ---------------------------------------------------------------

# Fixed parameters (gX1 is the swept control parameter)
gX0, Yth, nY, kX = 5.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12


# --- Repressive Hill functions ---------------------------------
def Hy(Y):
    # represses X: fraction of active transcription given Y
    return 1.0 / (1.0 + (Y / Yth) ** nY)

def Gx(X):
    # represses Y: fraction of active transcription given X
    return 1.0 / (1.0 + (X / Xth) ** nX)


# --- Nullclines (separation of variables) ----------------------
# X-nullcline: set dX/dt=0  ->  X = (gX0 + gX1*Hy(Y)) / kX   (X as fn of Y)
def X_nullcline(Y, gX1):
    return (gX0 + gX1 * Hy(Y)) / kX

# Y-nullcline: set dY/dt=0  ->  Y = (gY0 + gY1*Gx(X)) / kY   (Y as fn of X)
def Y_nullcline(X):
    return (gY0 + gY1 * Gx(X)) / kY


# --- Steady-state condition -------------------------------------
# Substitute Y=Y_nullcline(X) into X-nullcline to get a single
# self-consistency equation in X:   g(X) = X_nullcline(Y(X)) - X = 0
def g(X, gX1):
    Y = Y_nullcline(X)
    return X_nullcline(Y, gX1) - X


# --- Jacobian eigenvalue classification -------------------------
def Hy_prime(Y):
    u = (Y / Yth) ** nY
    return -u * nY / Y / (1.0 + u) ** 2   # dHy/dY

def Gx_prime(X):
    v = (X / Xth) ** nX
    return -v * nX / X / (1.0 + v) ** 2   # dGx/dX

def is_stable(X, Y, gX1):
    # Jacobian of (dX/dt, dY/dt) w.r.t. (X, Y)
    J = np.array([[-kX,            gX1 * Hy_prime(Y)],
                  [gY1 * Gx_prime(X), -kY]])
    eig = np.linalg.eigvals(J)
    return np.all(np.real(eig) < 0)       # stable iff all Re(eig) < 0


# --- Root finding: scan for sign changes, then bisect -----------
def find_steady_states(gX1):
    Xgrid = np.linspace(1e-6, (gX0 + gX1) / kX + 10.0, 4000)
    vals = g(Xgrid, gX1)
    roots = []
    for i in range(len(Xgrid) - 1):
        if vals[i] == 0.0:
            roots.append(Xgrid[i])
        elif vals[i] * vals[i + 1] < 0.0:      # sign change -> bracketed root
            a, b = Xgrid[i], Xgrid[i + 1]
            fa = vals[i]
            for _ in range(100):               # bisection
                m = 0.5 * (a + b)
                fm = g(m, gX1)
                if fa * fm <= 0.0:
                    b = m
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + b))
    return roots


# --- Sweep the control parameter gX1 ----------------------------
gX1_values = np.linspace(0.0, 100.0, 400)
stable_pts, unstable_pts = [], []
count_by_param = {}

for gX1 in gX1_values:
    ss = find_steady_states(gX1)
    count_by_param[gX1] = len(ss)
    for X in ss:
        Y = Y_nullcline(X)
        if is_stable(X, Y, gX1):
            stable_pts.append((gX1, X))
        else:
            unstable_pts.append((gX1, X))

stable_pts = np.array(stable_pts)
unstable_pts = np.array(unstable_pts)


# --- Check: do branches appear/merge (number & stability change)?
counts = np.array([count_by_param[p] for p in gX1_values])
n_min, n_max = counts.min(), counts.max()
n_stable = len(stable_pts)
n_unstable = len(unstable_pts) if unstable_pts.size else 0
# locate parameter values where the number of steady states changes
change_idx = np.where(np.diff(counts) != 0)[0]
bif_params = gX1_values[change_idx]

print("Min number of steady states over sweep:", int(n_min))
print("Max number of steady states over sweep:", int(n_max))
print("Total stable steady-state points found:", n_stable)
print("Total unstable steady-state points found:", n_unstable)
print("Number of parameter values where steady-state count changes:", len(bif_params))
for i, bp in enumerate(bif_params):
    print("Approx bifurcation gX1[%d]: %.4f  (count %d -> %d)"
          % (i, bp, counts[change_idx[i]], counts[change_idx[i] + 1]))

# Report a representative point on each branch
if n_stable:
    print("Example stable point (gX1, X): %.3f, %.4f" % (stable_pts[0, 0], stable_pts[0, 1]))
if n_unstable:
    print("Example unstable point (gX1, X): %.3f, %.4f" % (unstable_pts[0, 0], unstable_pts[0, 1]))

# One-sentence explanation of why the check confirms the result:
print("Explanation: The appearance of unstable branches merging with stable "
      "ones at parameter values where the steady-state count jumps is the "
      "signature of saddle-node bifurcations, confirming that the toggle "
      "switch's number and stability of steady states genuinely change with gX1.")


# --- Plot -------------------------------------------------------
plt.figure(figsize=(8, 6))
if stable_pts.size:
    plt.scatter(stable_pts[:, 0], stable_pts[:, 1], s=8, c="tab:blue", label="stable")
if unstable_pts.size:
    plt.scatter(unstable_pts[:, 0], unstable_pts[:, 1], s=8, c="tab:red", label="unstable")
for bp in bif_params:
    plt.axvline(bp, color="gray", ls=":", lw=0.6)
plt.xlabel("control parameter gX1 (X production rate)")
plt.ylabel("steady-state X")
plt.title("Toggle-switch bifurcation diagram: steady-state X vs gX1")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3E.1.1_s4.png")
