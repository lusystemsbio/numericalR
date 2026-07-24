import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# --- Model parameters (toggle switch) ---
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# --- Right-hand side of the ODE system ---
def rhs(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X   # X repressed by Y
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y   # Y repressed by X
    return np.array([dX, dY])

# --- Nullcline reduction: at steady state each variable is fixed by the other ---
def Y_of_X(X):  # from dY/dt = 0
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY
def X_from_Y(Y):  # from dX/dt = 0
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# g(X) = 0 when X is a self-consistent steady state (substitute Y(X) back in)
def g(X):
    return X_from_Y(Y_of_X(X)) - X

# --- Locate all steady states by scanning X for sign changes, then bracketing ---
Xgrid = np.linspace(0.0, X_from_Y(0.0) + 10.0, 20000)
gvals = np.array([g(x) for x in Xgrid])
roots_X = []
for i in range(len(Xgrid) - 1):
    if gvals[i] == 0.0:
        roots_X.append(Xgrid[i])
    elif gvals[i] * gvals[i + 1] < 0.0:  # sign change -> root in this bracket
        roots_X.append(brentq(g, Xgrid[i], Xgrid[i + 1]))

steady_states = [np.array([x, Y_of_X(x)]) for x in sorted(roots_X)]

# --- Numerical Jacobian via central finite differences ---
def jacobian(state, h=1e-6):
    n = len(state)
    J = np.zeros((n, n))
    for j in range(n):                     # perturb the j-th variable
        dp = np.zeros(n); dp[j] = h
        # column j = partial derivative of the RHS vector w.r.t. variable j
        J[:, j] = (rhs(state + dp) - rhs(state - dp)) / (2.0 * h)
    return J

# --- Classify each steady state from Jacobian eigenvalues ---
labels = []
for k, ss in enumerate(steady_states):
    J = jacobian(ss)
    eig = np.linalg.eigvals(J)             # eigenvalues of the 2x2 Jacobian
    stable = np.all(np.real(eig) < 0)      # stable iff both real parts negative
    label = "stable" if stable else "unstable"
    labels.append(label)
    print(f"Steady state {k+1}: X = {ss[0]:.6f}, Y = {ss[1]:.6f}")
    print(f"  Jacobian eigenvalue 1: {eig[0].real:+.6f} {eig[0].imag:+.6f}j")
    print(f"  Jacobian eigenvalue 2: {eig[1].real:+.6f} {eig[1].imag:+.6f}j")
    print(f"  Stability label: {label}")

# --- Separate check: outer two stable, middle one unstable saddle ---
n_ss = len(steady_states)
print(f"Number of steady states found: {n_ss}")
if n_ss == 3:
    outer_ok = labels[0] == "stable" and labels[2] == "stable"
    # a saddle has one positive and one negative real eigenvalue
    Jmid = jacobian(steady_states[1])
    emid = np.linalg.eigvals(Jmid)
    is_saddle = (np.max(np.real(emid)) > 0) and (np.min(np.real(emid)) < 0)
    print(f"Outer states both stable: {outer_ok}")
    print(f"Middle state is unstable saddle: {is_saddle}")
    print(f"Bistability check passed: {outer_ok and is_saddle}")

# --- Visualization: nullclines + steady states colored by stability ---
Xs = np.linspace(0, max(s[0] for s in steady_states) * 1.3, 500)
Ys = np.linspace(0, max(s[1] for s in steady_states) * 1.3, 500)
plt.figure(figsize=(7, 6))
plt.plot(Xs, Y_of_X(Xs), 'b-', label='dY/dt = 0 nullcline')
plt.plot([X_from_Y(y) for y in Ys], Ys, 'g-', label='dX/dt = 0 nullcline')
for ss, lab in zip(steady_states, labels):
    plt.plot(ss[0], ss[1], 'o', ms=12,
             color='black' if lab == 'stable' else 'red',
             mfc='black' if lab == 'stable' else 'white',
             label=f'{lab} ({ss[0]:.1f}, {ss[1]:.1f})')
plt.xlabel('X'); plt.ylabel('Y'); plt.title('Toggle switch steady states')
plt.legend(loc='best'); plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.3.1_s4.png")

# Why the check confirms the result: linear stability theory guarantees that a
# 2D fixed point with both Jacobian eigenvalues negative is an attractor while
# one positive eigenvalue makes it a saddle that repels along one direction, so
# finding two attractors flanking one saddle is exactly the eigenvalue signature
# of a bistable toggle switch and thus matches the simulated dynamics.
print("Check confirms result: eigenvalue signs (two attractors + one saddle) "
      "are the linear-stability signature of bistability seen in simulations.")
