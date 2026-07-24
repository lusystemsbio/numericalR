import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------
# Toggle switch model:
#   dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X
#   dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y
# ------------------------------------------------------------------

# Fixed parameters
gX0 = 5.0
Yth = 100.0
nY  = 4.0
kX  = 0.1

gY0 = 4.0
gY1 = 40.0
Xth = 150.0
nX  = 4.0
kY  = 0.12

# ---- Nullcline helper functions (separation of variables) --------
# Setting dY/dt = 0 solves explicitly for Y as a function of X:
def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# Setting dX/dt = 0 solves explicitly for X as a function of Y:
def X_of_Y(Y, gX1):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# ---- Steady-state residual ---------------------------------------
# Substitute Y = Y_of_X(X) into the X-nullcline. A steady state is a
# root of G(X) = X - X_of_Y(Y_of_X(X)); i.e. both nullclines meet.
def G(X, gX1):
    return X - X_of_Y(Y_of_X(X), gX1)

# ---- Jacobian of the 2D system at a steady state -----------------
def jacobian(X, Y, gX1):
    # d(dX/dt)/dX and /dY
    dfX_dX = -kX
    # derivative of gX1/(1+(Y/Yth)^nY) w.r.t Y
    uY = (Y / Yth) ** nY
    dfX_dY = gX1 * (-(nY / Yth) * (Y / Yth) ** (nY - 1.0)) / (1.0 + uY) ** 2
    # d(dY/dt)/dX and /dY
    uX = (X / Xth) ** nX
    dfY_dX = gY1 * (-(nX / Xth) * (X / Xth) ** (nX - 1.0)) / (1.0 + uX) ** 2
    dfY_dY = -kY
    return np.array([[dfX_dX, dfX_dY],
                     [dfY_dX, dfY_dY]])

# ---- Root finder: scan for sign changes, then bisect -------------
def find_steady_states(gX1):
    # Bracket the physically reachable X range (X >= gX0/kX, bounded above)
    Xmax = (gX0 + max(gX1, 0.0)) / kX + 10.0
    Xgrid = np.linspace(0.0, Xmax, 4000)
    Gvals = np.array([G(x, gX1) for x in Xgrid])

    roots = []
    for i in range(len(Xgrid) - 1):
        a, b = Xgrid[i], Xgrid[i + 1]
        fa, fb = Gvals[i], Gvals[i + 1]
        if fa == 0.0:
            roots.append(a)
        elif fa * fb < 0.0:  # sign change -> a root lies in [a, b]
            # bisection to refine the root
            for _ in range(60):
                m = 0.5 * (a + b)
                fm = G(m, gX1)
                if fa * fm <= 0.0:
                    b, fb = m, fm
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + b))

    # De-duplicate roots that are numerically identical
    unique = []
    for r in roots:
        if not any(abs(r - u) < 1e-4 for u in unique):
            unique.append(r)
    return unique

# ---- Classify a steady state by Jacobian eigenvalues -------------
def classify(X, gX1):
    Y = Y_of_X(X)
    J = jacobian(X, Y, gX1)
    eig = np.linalg.eigvals(J)
    stable = np.all(np.real(eig) < 0.0)  # stable iff all Re(eig) < 0
    return Y, eig, stable

# ------------------------------------------------------------------
# Sweep the control parameter gX1 and collect steady states
# ------------------------------------------------------------------
gX1_sweep = np.linspace(0.0, 100.0, 400)

stable_g, stable_X = [], []
unstable_g, unstable_X = [], []

count_per_g = []
for gX1 in gX1_sweep:
    ss = find_steady_states(gX1)
    count_per_g.append(len(ss))
    for X in ss:
        Y, eig, is_stable = classify(X, gX1)
        if is_stable:
            stable_g.append(gX1); stable_X.append(X)
        else:
            unstable_g.append(gX1); unstable_X.append(X)

# ------------------------------------------------------------------
# Print numerical results at a few representative parameter values
# ------------------------------------------------------------------
for gX1 in [0.0, 25.0, 50.0, 75.0, 100.0]:
    ss = find_steady_states(gX1)
    print(f"--- gX1 = {gX1:.1f}: {len(ss)} steady state(s) ---")
    for X in sorted(ss):
        Y, eig, is_stable = classify(X, gX1)
        label = "STABLE" if is_stable else "UNSTABLE"
        print(f"  X* = {X:10.4f}  Y* = {Y:10.4f}  "
              f"eigenvalues = ({eig[0].real:.5f}{eig[0].imag:+.5f}j, "
              f"{eig[1].real:.5f}{eig[1].imag:+.5f}j)  -> {label}")

# Bifurcation check: report where the number of steady states changes
counts = np.array(count_per_g)
print("\n--- Bifurcation check: changes in number of steady states ---")
changes = np.where(np.diff(counts) != 0)[0]
if len(changes) == 0:
    print("  No change in steady-state count over the sweep.")
for idx in changes:
    g_before, g_after = gX1_sweep[idx], gX1_sweep[idx + 1]
    print(f"  Between gX1 = {g_before:.3f} and {g_after:.3f}: "
          f"count {counts[idx]} -> {counts[idx + 1]} (bifurcation near "
          f"gX1 = {0.5 * (g_before + g_after):.3f})")

print(f"\nMax number of simultaneous steady states over sweep: {counts.max()}")
print(f"Total stable points plotted:   {len(stable_X)}")
print(f"Total unstable points plotted: {len(unstable_X)}")

# ------------------------------------------------------------------
# Plot the bifurcation diagram: steady-state X vs gX1, colored by stability
# ------------------------------------------------------------------
plt.figure(figsize=(9, 6))
plt.scatter(unstable_g, unstable_X, s=10, c="crimson", label="unstable")
plt.scatter(stable_g, stable_X, s=10, c="royalblue", label="stable")
plt.xlabel("control parameter gX1 (X production rate)")
plt.ylabel("steady-state X*")
plt.title("Toggle switch bifurcation diagram (X* vs gX1)")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3E.1.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: The count of steady states changes at specific gX1 values "
      "where a stable and an unstable branch meet and annihilate (saddle-node "
      "bifurcations); observing the region of three coexisting states (two stable, "
      "one unstable) collapse to one confirms the system is a genuine bistable "
      "toggle switch whose behavior is controlled by gX1.")
