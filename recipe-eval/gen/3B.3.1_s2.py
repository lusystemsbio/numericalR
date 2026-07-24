import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# ---- Model parameters (toggle switch: X and Y repress each other) ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ---- Right-hand side of the ODE system ----
def f(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X   # X production - degradation
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y   # Y production - degradation
    return np.array([dX, dY])

# ---- Find steady states: solve f=0 from many initial guesses, keep unique roots ----
guesses = [(x, y) for x in np.linspace(0, 600, 12) for y in np.linspace(0, 600, 12)]
roots = []
for g in guesses:
    sol, info, ier, msg = fsolve(f, g, full_output=True)
    if ier == 1 and np.max(np.abs(f(sol))) < 1e-9:            # converged to a real root
        if not any(np.allclose(sol, r, atol=1e-4) for r in roots):  # keep only new ones
            roots.append(sol)
roots = sorted(roots, key=lambda r: r[0])   # order by X to identify low-X / mid / high-X

# ---- Numerical Jacobian via central finite differences ----
def jacobian(state, h=1e-6):
    J = np.zeros((2, 2))
    for j in range(2):                       # perturb each variable j
        dp = np.array(state, dtype=float); dp[j] += h
        dm = np.array(state, dtype=float); dm[j] -= h
        J[:, j] = (f(dp) - f(dm)) / (2 * h)  # column j = d f / d state_j
    return J

# ---- Classify each steady state from the eigenvalues of the Jacobian ----
print("Found %d steady states" % len(roots))
labels = []
for i, ss in enumerate(roots):
    J = jacobian(ss)
    eig = np.linalg.eigvals(J)               # eigenvalues of the 2x2 linearization
    stable = np.all(np.real(eig) < 0)        # stable iff BOTH real parts are negative
    saddle = (np.real(eig)[0] * np.real(eig)[1]) < 0  # real parts of opposite sign
    label = "STABLE" if stable else ("UNSTABLE (saddle)" if saddle else "UNSTABLE")
    labels.append(label)
    print("----------------------------------------")
    print("Steady state %d: X = %.6f, Y = %.6f" % (i + 1, ss[0], ss[1]))
    print("  eigenvalue 1 = %s" % np.format_float_positional(eig[0].real, precision=6) +
          ("+%.6gj" % eig[0].imag if eig[0].imag else ""))
    print("  eigenvalue 2 = %s" % np.format_float_positional(eig[1].real, precision=6) +
          ("+%.6gj" % eig[1].imag if eig[1].imag else ""))
    print("  Re(eig): %.6f, %.6f" % (eig[0].real, eig[1].real))
    print("  Stability label: %s" % label)

# ---- Separate check: outer states stable, middle is an unstable saddle ----
print("========================================")
outer_ok = (labels[0] == "STABLE") and (labels[-1] == "STABLE")
mid_ok = ("saddle" in labels[len(labels) // 2].lower())
print("Outer states both STABLE:            %s" % outer_ok)
print("Middle state UNSTABLE saddle:        %s" % mid_ok)
print("Matches expected bistable structure: %s" % (outer_ok and mid_ok))
# One-sentence explanation:
print("Explanation: two stable fixed points flanking a single saddle is exactly the "
      "phase-portrait signature of bistability, so this eigenvalue-based classification "
      "reproduces what the toggle-switch simulations show.")

# ---- Visualize: nullclines, steady states colored by stability ----
Xg = np.linspace(0, 600, 400)
Yg = np.linspace(0, 600, 400)
# X-nullcline (dX=0): X as function of Y ; Y-nullcline (dY=0): Y as function of X
X_null = (gX0 + gX1 / (1.0 + (Yg / Yth) ** nY)) / kX   # gives X for each Y
Y_null = (gY0 + gY1 / (1.0 + (Xg / Xth) ** nX)) / kY   # gives Y for each X
plt.figure(figsize=(6, 5))
plt.plot(X_null, Yg, label="dX/dt = 0 nullcline")
plt.plot(Xg, Y_null, label="dY/dt = 0 nullcline")
for ss, lab in zip(roots, labels):
    color = "green" if lab == "STABLE" else "red"
    plt.plot(ss[0], ss[1], "o", color=color, markersize=10,
             markeredgecolor="black", zorder=5)
    plt.annotate(lab, (ss[0], ss[1]), textcoords="offset points", xytext=(8, 8))
plt.xlabel("X"); plt.ylabel("Y"); plt.title("Toggle switch steady states")
plt.legend(); plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.3.1_s2.png")
print("Saved figure.")
