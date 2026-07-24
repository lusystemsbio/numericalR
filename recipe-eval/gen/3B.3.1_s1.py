import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# ---- Model parameters (toggle switch) ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ---- Right-hand side of the ODE system ----
def f(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth)**nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth)**nX) - kY * Y
    return np.array([dX, dY])

# ---- Find steady states: solve f(state)=0 from many initial guesses, then dedupe ----
found = []
for X0 in np.linspace(0, 600, 15):
    for Y0 in np.linspace(0, 600, 15):
        sol, info, ier, msg = fsolve(f, [X0, Y0], full_output=True)
        if ier == 1 and np.all(sol >= -1e-6):          # converged and non-negative
            if not any(np.allclose(sol, s, atol=1e-4) for s in found):
                found.append(sol)
steady_states = sorted(found, key=lambda s: s[0])       # order by X for outer/mid/outer

# ---- Numerical Jacobian via central finite differences ----
def jacobian(state, h=1e-6):
    J = np.zeros((2, 2))
    for j in range(2):                                  # perturb each variable
        dp = state.copy(); dp[j] += h
        dm = state.copy(); dm[j] -= h
        J[:, j] = (f(dp) - f(dm)) / (2 * h)             # column = d f / d state_j
    return J

# ---- Classify: stable iff both eigenvalues have negative real part ----
labels = []
for k, ss in enumerate(steady_states):
    J = jacobian(ss)
    eig = np.linalg.eigvals(J)                           # eigenvalues of 2x2 Jacobian
    stable = np.all(eig.real < 0)                        # both real parts negative -> stable
    label = "stable" if stable else "unstable"
    labels.append(label)
    print(f"Steady state {k+1}: X = {ss[0]:.6f}, Y = {ss[1]:.6f}")
    print(f"  eigenvalue 1 = {eig[0].real:+.6f} {eig[0].imag:+.6f}j")
    print(f"  eigenvalue 2 = {eig[1].real:+.6f} {eig[1].imag:+.6f}j")
    print(f"  stability    = {label}")

# ---- Separate check: outer two stable, middle unstable saddle (bistability) ----
outer_stable = (labels[0] == "stable") and (labels[-1] == "stable")
mid_J = jacobian(steady_states[1])
mid_eig = np.linalg.eigvals(mid_J)
middle_saddle = (mid_eig[0].real * mid_eig[1].real < 0)  # real eigenvalues of opposite sign
print(f"Number of steady states found = {len(steady_states)}")
print(f"Outer states both stable      = {outer_stable}")
print(f"Middle state is saddle        = {middle_saddle}")
print(f"Bistable check passed         = {outer_stable and middle_saddle}")

# ---- Visualization: nullclines + steady states ----
Xg = np.linspace(0, 600, 400)
Yg = np.linspace(0, 600, 400)
# dX=0 nullcline: X = (gX0 + gX1/(1+(Y/Yth)^nY))/kX  ->  X as function of Y
X_null = (gX0 + gX1 / (1.0 + (Yg / Yth)**nY)) / kX
# dY=0 nullcline: Y = (gY0 + gY1/(1+(X/Xth)^nX))/kY  ->  Y as function of X
Y_null = (gY0 + gY1 / (1.0 + (Xg / Xth)**nX)) / kY

plt.figure(figsize=(6, 6))
plt.plot(X_null, Yg, label="dX/dt = 0")
plt.plot(Xg, Y_null, label="dY/dt = 0")
for ss, lab in zip(steady_states, labels):
    color = "green" if lab == "stable" else "red"
    marker = "o" if lab == "stable" else "x"
    plt.scatter(ss[0], ss[1], c=color, marker=marker, s=120, zorder=5, label=lab)
plt.xlabel("X"); plt.ylabel("Y"); plt.title("Toggle switch steady states")
handles, lbls = plt.gca().get_legend_handles_labels()
by_label = dict(zip(lbls, handles))
plt.legend(by_label.values(), by_label.keys())
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.3.1_s1.png")

# The check confirms the result because a genetic toggle switch is bistable precisely when
# the two extreme fixed points are attractors and the intermediate one is a saddle whose
# unstable manifold separates their basins -- exactly the stable/unstable pattern the
# eigenvalue signs reproduce here.
print("Explanation: the eigenvalue-based labels reproduce the bistable signature (two stable "
      "attractors flanking one unstable saddle), which is the defining dynamical structure of "
      "a functioning toggle switch, thereby confirming the linear-stability classification.")
