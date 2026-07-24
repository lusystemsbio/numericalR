import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Toggle-switch parameters ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ---- Closed-form nullclines by separation of variables ----
# X-nullcline: set dX/dt = 0 and solve for X (X is isolated linearly):
#   gX0 + gX1/(1+(Y/Yth)^nY) - kX*X = 0  ->  X = [gX0 + gX1/(1+(Y/Yth)^nY)] / kX
def X_of_Y(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# Y-nullcline: set dY/dt = 0 and solve for Y (Y is isolated linearly):
#   gY0 + gY1/(1+(X/Xth)^nX) - kY*Y = 0  ->  Y = [gY0 + gY1/(1+(X/Xth)^nX)] / kY
def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ---- Sweep the "other" variable to trace each curve ----
Y_sweep = np.linspace(0.0, 600.0, 2000)   # for X-nullcline: X as function of Y
X_null_X = X_of_Y(Y_sweep)                 # X coordinates of the X-nullcline

X_sweep = np.linspace(0.0, 600.0, 2000)   # for Y-nullcline: Y as function of X
Y_null_Y = Y_of_X(X_sweep)                 # Y coordinates of the Y-nullcline

# ---- Find intersections (steady states) ----
# A steady state satisfies both simultaneously. Compose the two closed forms:
#   X* must satisfy X* = X_of_Y( Y_of_X(X*) ).  Define g(X) = X_of_Y(Y_of_X(X)) - X
# and locate sign changes, then refine each bracket by bisection.
def g(X):
    return X_of_Y(Y_of_X(X)) - X

Xgrid = np.linspace(1e-6, 600.0, 20000)   # scan for sign changes
gvals = g(Xgrid)
brackets = []
for i in range(len(Xgrid) - 1):
    if gvals[i] == 0.0:
        brackets.append((Xgrid[i], Xgrid[i]))
    elif gvals[i] * gvals[i + 1] < 0.0:
        brackets.append((Xgrid[i], Xgrid[i + 1]))

steady_states = []
for a, b in brackets:
    # simple bisection to refine the root of g
    for _ in range(80):
        m = 0.5 * (a + b)
        if g(a) * g(m) <= 0.0:
            b = m
        else:
            a = m
    Xs = 0.5 * (a + b)
    Ys = Y_of_X(Xs)
    steady_states.append((Xs, Ys))

# ---- Report numerical results ----
print(f"Number of nullcline crossings (steady states) found: {len(steady_states)}")
for idx, (Xs, Ys) in enumerate(steady_states, start=1):
    # residuals of the original ODE right-hand sides, should be ~0
    fX = gX0 + gX1 / (1.0 + (Ys / Yth) ** nY) - kX * Xs
    fY = gY0 + gY1 / (1.0 + (Xs / Xth) ** nX) - kY * Ys
    print(f"Steady state {idx}: X = {Xs:.6f}, Y = {Ys:.6f}  (residual fX = {fX:.3e}, fY = {fY:.3e})")

# ---- Stability check via linearization (Jacobian eigenvalues) ----
def jacobian(X, Y):
    # d(fX)/dX = -kX ; d(fX)/dY = gX1 * d/dY[1/(1+(Y/Yth)^nY)]
    dfX_dX = -kX
    dHY = -nY * (Y / Yth) ** (nY - 1) / Yth / (1.0 + (Y / Yth) ** nY) ** 2
    dfX_dY = gX1 * dHY
    # d(fY)/dY = -kY ; d(fY)/dX = gY1 * d/dX[1/(1+(X/Xth)^nX)]
    dfY_dY = -kY
    dHX = -nX * (X / Xth) ** (nX - 1) / Xth / (1.0 + (X / Xth) ** nX) ** 2
    dfY_dX = gY1 * dHX
    return np.array([[dfX_dX, dfX_dY], [dfY_dX, dfY_dY]])

n_stable, n_unstable = 0, 0
for idx, (Xs, Ys) in enumerate(steady_states, start=1):
    eig = np.linalg.eigvals(jacobian(Xs, Ys))
    stable = np.all(eig.real < 0.0)
    if stable:
        n_stable += 1
    else:
        n_unstable += 1
    print(f"Steady state {idx}: eigenvalues = {eig[0]:.4f}, {eig[1]:.4f}  -> {'STABLE' if stable else 'UNSTABLE'}")

print(f"Stable steady states:   {n_stable}")
print(f"Unstable steady states: {n_unstable}")
print(f"Three-crossing check passes: {len(steady_states) == 3 and n_stable == 2 and n_unstable == 1}")

# ---- Phase-plane plot of both nullclines and their crossings ----
plt.figure(figsize=(7, 6))
plt.plot(X_null_X, Y_sweep, 'b-', label="X-nullcline (dX/dt = 0)")
plt.plot(X_sweep, Y_null_Y, 'r-', label="Y-nullcline (dY/dt = 0)")
for idx, (Xs, Ys) in enumerate(steady_states, start=1):
    eig = np.linalg.eigvals(jacobian(Xs, Ys))
    stable = np.all(eig.real < 0.0)
    plt.plot(Xs, Ys, 'ko' if stable else 'ks',
             markersize=10, markerfacecolor=('k' if stable else 'white'),
             label=("stable steady state" if (stable and idx == 1) else
                    ("unstable steady state" if not stable else None)))
    plt.annotate(f"({Xs:.0f}, {Ys:.0f})", (Xs, Ys),
                 textcoords="offset points", xytext=(8, 8))
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Toggle switch nullclines and steady states")
plt.legend()
plt.xlim(0, 600)
plt.ylim(0, 600)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.4.1_s2.png")

# Explanation: the two nullclines are the loci where each variable stops changing,
# so every point where they cross is a steady state; counting three crossings that
# split into two stable and one unstable node confirms the bistable toggle-switch
# behavior expected from direct time-domain simulation.
print("Why the check confirms it: intersections of the nullclines are exactly the")
print("points where both dX/dt=0 and dY/dt=0, so three crossings (two stable, one")
print("unstable by Jacobian eigenvalues) reproduce the bistable steady states a")
print("time-domain simulation would settle into.")
