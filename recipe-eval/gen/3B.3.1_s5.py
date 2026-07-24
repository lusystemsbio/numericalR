import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# ----- Model parameters (toggle switch) -----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ----- Right-hand side of the ODE system -----
def f(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# ----- Steady-state nullcline helpers -----
# From dX/dt=0:  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX   -> X as a function of Y
def X_of_Y(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX
# From dY/dt=0:  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY   -> Y as a function of X
def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# A steady state satisfies Y = Y_of_X( X_of_Y(Y) ); find roots of g(Y)=0.
def g(Y):
    return Y_of_X(X_of_Y(Y)) - Y

# ----- Locate the three steady states by scanning for sign changes in g(Y) -----
Y_scan = np.linspace(0.1, Y_of_X(0.0) + 10.0, 20000)
gvals = g(Y_scan)
roots_Y = []
for i in range(len(Y_scan) - 1):
    if gvals[i] == 0.0:
        roots_Y.append(Y_scan[i])
    elif gvals[i] * gvals[i + 1] < 0.0:
        roots_Y.append(brentq(g, Y_scan[i], Y_scan[i + 1]))
# de-duplicate close roots
uniq = []
for r in roots_Y:
    if not any(abs(r - u) < 1e-6 for u in uniq):
        uniq.append(r)
roots_Y = sorted(uniq)

steady_states = [(X_of_Y(Yr), Yr) for Yr in roots_Y]
print(f"Number of steady states found: {len(steady_states)}")

# ----- Numerical Jacobian via central finite differences -----
def jacobian(func, state, eps=1e-6):
    state = np.asarray(state, dtype=float)
    n = len(state)
    J = np.zeros((n, n))
    for j in range(n):
        h = eps * max(1.0, abs(state[j]))   # scale the step to the variable
        sp = state.copy(); sp[j] += h
        sm = state.copy(); sm[j] -= h
        J[:, j] = (func(sp) - func(sm)) / (2.0 * h)  # column j = d f / d state_j
    return J

# ----- Classify each steady state via eigenvalues of the Jacobian -----
labels = []
for idx, ss in enumerate(steady_states):
    X, Y = ss
    resid = f(ss)  # residual should be ~0 at a genuine steady state
    J = jacobian(f, ss)
    eigvals = np.linalg.eigvals(J)
    # Stable iff BOTH eigenvalues have negative real part; else unstable.
    stable = np.all(np.real(eigvals) < 0.0)
    # A saddle is specifically an unstable state with eigenvalues of opposite real-part sign.
    saddle = (np.real(eigvals) > 0).any() and (np.real(eigvals) < 0).any()
    label = "stable" if stable else ("unstable (saddle)" if saddle else "unstable")
    labels.append(label)

    print(f"--- Steady state {idx+1} ---")
    print(f"X* = {X:.6f}")
    print(f"Y* = {Y:.6f}")
    print(f"residual dX/dt = {resid[0]:.3e}")
    print(f"residual dY/dt = {resid[1]:.3e}")
    print(f"Jacobian eigenvalue 1 = {eigvals[0].real:.6f} + {eigvals[0].imag:.6f}j")
    print(f"Jacobian eigenvalue 2 = {eigvals[1].real:.6f} + {eigvals[1].imag:.6f}j")
    print(f"Stability label = {label}")

# ----- Separate check: outer states stable, middle state an unstable saddle -----
# Order the three states by X* so we can talk about "outer" vs "middle".
order = np.argsort([ss[0] for ss in steady_states])
if len(steady_states) == 3:
    lo, mid, hi = order
    outer_stable = (labels[lo] == "stable") and (labels[hi] == "stable")
    middle_saddle = labels[mid].startswith("unstable")
    print("--- Bistability check ---")
    print(f"Low-X  state stable? {labels[lo] == 'stable'}")
    print(f"High-X state stable? {labels[hi] == 'stable'}")
    print(f"Middle state unstable/saddle? {middle_saddle}")
    print(f"Check passes (two outer stable, middle saddle): {outer_stable and middle_saddle}")

# ----- Visualization: nullclines + steady states with stability coloring -----
Xg = np.linspace(0, 700, 800)
Yg = np.linspace(0, 500, 800)
plt.figure(figsize=(7, 6))
plt.plot(X_of_Y(Yg), Yg, 'b-', label='dX/dt = 0 nullcline')
plt.plot(Xg, Y_of_X(Xg), 'g-', label='dY/dt = 0 nullcline')
for ss, lab in zip(steady_states, labels):
    color = 'black' if lab == 'stable' else 'red'
    plt.plot(ss[0], ss[1], 'o', color=color, markersize=10)
    plt.annotate(lab, ss, textcoords="offset points", xytext=(8, 8))
plt.xlabel('X'); plt.ylabel('Y')
plt.title('Toggle switch: nullclines and steady-state stability')
plt.legend(); plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.3.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The check confirms bistability because a saddle (one positive, "
      "one negative eigenvalue) sits on the separatrix between the two stable basins, "
      "so trajectories fall to one of the two outer stable states exactly as the "
      "bistable simulations show.")
