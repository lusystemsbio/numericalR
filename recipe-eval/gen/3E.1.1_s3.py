import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle switch model (genes X and Y mutually repressing):
#   dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X
#   dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y
# ----------------------------------------------------------------------

# Fixed parameters
gX0 = 5.0;  Yth = 100.0; nY = 4.0; kX = 0.10
gY0 = 4.0;  gY1 = 40.0;  Xth = 150.0; nX = 4.0; kY = 0.12

# --- Nullclines via "separation of variables" -------------------------
# Setting each time derivative to zero and solving for that variable:
#   X-nullcline: X = ( gX0 + gX1/(1+(Y/Yth)^nY) ) / kX      -> X as function of Y
#   Y-nullcline: Y = ( gY0 + gY1/(1+(X/Xth)^nX) ) / kY      -> Y as function of X
def X_of_Y(Y, gX1):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# A steady state satisfies X = X_of_Y(Y_of_X(X)); define residual g(X).
def g(X, gX1):
    return X_of_Y(Y_of_X(X), gX1) - X

# --- Root finding: scan for sign changes, then bisect -----------------
def bisect(f, a, b, tol=1e-10, itmax=200):
    fa, fb = f(a), f(b)
    for _ in range(itmax):
        m = 0.5 * (a + b)
        fm = f(m)
        if abs(fm) < tol or 0.5 * (b - a) < tol:
            return m
        if (fa > 0) != (fm > 0):
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)

def find_steady_states(gX1):
    # X ranges from basal (gX0/kX) up to full production ((gX0+gX1)/kX)
    Xgrid = np.linspace(1e-6, (gX0 + gX1) / kX + 10.0, 4000)
    vals = np.array([g(x, gX1) for x in Xgrid])
    roots = []
    for i in range(len(Xgrid) - 1):
        if vals[i] == 0.0:
            roots.append(Xgrid[i])
        elif (vals[i] > 0) != (vals[i + 1] > 0):  # sign change -> a root between
            roots.append(bisect(lambda x: g(x, gX1), Xgrid[i], Xgrid[i + 1]))
    return roots

# --- Jacobian and stability classification ----------------------------
def jacobian(X, Y, gX1):
    # partial derivatives of the Hill terms
    # d/dY [ 1/(1+(Y/Yth)^nY) ] = -nY*(Y/Yth)^nY / (Y * (1+(Y/Yth)^nY)^2)
    uY = (Y / Yth) ** nY
    dFx_dY = gX1 * (-nY * uY / (Y * (1.0 + uY) ** 2))
    uX = (X / Xth) ** nX
    dFy_dX = gY1 * (-nX * uX / (X * (1.0 + uX) ** 2))
    return np.array([[-kX, dFx_dY],
                     [dFy_dX, -kY]])

def is_stable(X, Y, gX1):
    ev = np.linalg.eigvals(jacobian(X, Y, gX1))
    return np.all(ev.real < 0)  # stable iff all eigenvalues have negative real part

# --- Sweep the control parameter gX1 ----------------------------------
gX1_values = np.linspace(0.0, 100.0, 400)
stab_g, stab_X = [], []
unst_g, unst_X = [], []
counts = []
for gX1 in gX1_values:
    ss = find_steady_states(gX1)
    counts.append(len(ss))
    for X in ss:
        Y = Y_of_X(X)
        if is_stable(X, Y, gX1):
            stab_g.append(gX1); stab_X.append(X)
        else:
            unst_g.append(gX1); unst_X.append(X)

# --- Print numerical results ------------------------------------------
print(f"Parameter sweep: gX1 from {gX1_values[0]:.1f} to {gX1_values[-1]:.1f} ({len(gX1_values)} points)")
print(f"Total stable steady-state points found: {len(stab_X)}")
print(f"Total unstable steady-state points found: {len(unst_X)}")

counts = np.array(counts)
print(f"Distinct steady-state counts observed across sweep: {sorted(set(counts.tolist()))}")
print(f"Max number of coexisting steady states: {counts.max()}")
print(f"Min number of coexisting steady states: {counts.min()}")

# Detect where the number of steady states changes (bifurcation points)
print("Bifurcation points (gX1 where the steady-state count changes):")
for i in range(1, len(gX1_values)):
    if counts[i] != counts[i - 1]:
        gb = 0.5 * (gX1_values[i] + gX1_values[i - 1])
        print(f"  gX1 ~= {gb:.3f} : count {counts[i-1]} -> {counts[i]}")

# Report a detailed steady-state table at a few sample gX1 values
for gX1 in [0.0, 25.0, 50.0, 75.0, 100.0]:
    ss = find_steady_states(gX1)
    labels = []
    for X in ss:
        Y = Y_of_X(X)
        labels.append(f"X={X:.3f}({'stable' if is_stable(X, Y, gX1) else 'unstable'})")
    print(f"gX1={gX1:6.1f}: {len(ss)} steady state(s): " + ", ".join(labels))

# --- Plot: steady-state X vs control parameter, colored by stability ---
plt.figure(figsize=(8, 6))
plt.scatter(unst_g, unst_X, s=10, c="red", label="unstable")
plt.scatter(stab_g, stab_X, s=10, c="blue", label="stable")
plt.xlabel("control parameter gX1 (X production rate)")
plt.ylabel("steady-state X")
plt.title("Toggle switch bifurcation diagram (colored by stability)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3E.1.1_s3.png")

# --- One-sentence explanation of the check ----------------------------
print("Explanation: Seeing stable and unstable branches that appear/merge at "
      "parameter values where the steady-state count changes confirms the result, "
      "because a saddle-node bifurcation is precisely the collision of a stable and "
      "an unstable fixed point, so the branch structure directly reflects the "
      "changing number and stability of steady states.")
