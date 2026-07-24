import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import fsolve

# ------------------------------------------------------------------
# Model parameters (toggle switch: genes X and Y mutually repress)
# ------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Right-hand side of the ODE system: f(s) = [dX/dt, dY/dt]
def f(s):
    X, Y = s
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# ------------------------------------------------------------------
# Numerical Jacobian by finite differences (central differences),
# implemented explicitly rather than via an autodiff/library routine.
# J[i,j] = d f_i / d s_j
# ------------------------------------------------------------------
def jacobian(s, h=1e-6):
    s = np.asarray(s, dtype=float)
    n = s.size
    J = np.zeros((n, n))
    for j in range(n):                      # perturb each variable
        step = np.zeros(n)
        step[j] = h
        f_plus = f(s + step)                # forward perturbation
        f_minus = f(s - step)               # backward perturbation
        J[:, j] = (f_plus - f_minus) / (2.0 * h)   # central difference
    return J

# ------------------------------------------------------------------
# Stability test: linearize about the steady state, take eigenvalues
# of the 2x2 Jacobian; stable iff every eigenvalue has Re < 0.
# ------------------------------------------------------------------
def classify(s):
    J = jacobian(s)
    eigvals = np.linalg.eigvals(J)          # eigenvalues of linearization
    stable = np.all(eigvals.real < 0.0)     # both real parts negative -> stable
    return J, eigvals, stable

# ------------------------------------------------------------------
# Find the three steady states by solving f(s) = 0 from several
# initial guesses spanning the (X, Y) plane, then de-duplicate.
# ------------------------------------------------------------------
guesses = [(50, 400), (250, 40), (150, 200), (10, 500), (500, 10),
           (100, 100), (200, 100), (100, 300)]
found = []
for g in guesses:
    sol, info, ier, msg = fsolve(f, g, full_output=True)
    if ier == 1 and np.allclose(f(sol), 0.0, atol=1e-6) and np.all(sol > 0):
        if not any(np.allclose(sol, s, atol=1e-3) for s in found):
            found.append(sol)

# Sort by X coordinate so "outer" and "middle" are well defined
found.sort(key=lambda s: s[0])

print("Number of steady states found:", len(found))
print()

results = []
for i, s in enumerate(found):
    J, eig, stable = classify(s)
    label = "STABLE" if stable else "UNSTABLE"
    results.append((s, eig, stable))
    print(f"Steady state {i+1}: X = {s[0]:.6f}, Y = {s[1]:.6f}")
    print(f"  Jacobian =\n{J}")
    print(f"  Eigenvalue 1: {eig[0].real:+.6f} {eig[0].imag:+.6f}j")
    print(f"  Eigenvalue 2: {eig[1].real:+.6f} {eig[1].imag:+.6f}j")
    print(f"  Stability label: {label}")
    # A saddle has one positive and one negative real eigenvalue
    if not stable and np.any(eig.real > 0) and np.any(eig.real < 0):
        print("  (unstable SADDLE: eigenvalues of opposite sign)")
    print()

# ------------------------------------------------------------------
# Separate check: outer two states stable, middle one unstable saddle
# ------------------------------------------------------------------
if len(results) == 3:
    outer_lo_stable = results[0][2]
    middle_stable = results[1][2]
    outer_hi_stable = results[2][2]
    middle_saddle = (not middle_stable
                     and np.any(results[1][1].real > 0)
                     and np.any(results[1][1].real < 0))
    print("CHECK - lower outer state stable :", outer_lo_stable)
    print("CHECK - upper outer state stable :", outer_hi_stable)
    print("CHECK - middle state is saddle   :", middle_saddle)
    check_pass = outer_lo_stable and outer_hi_stable and middle_saddle
    print("CHECK - matches bistable expectation:", check_pass)
else:
    print("CHECK - expected 3 steady states, cannot run bistability check")

# ------------------------------------------------------------------
# Visualization: nullclines, steady states, and stability coloring
# ------------------------------------------------------------------
X = np.linspace(0, 600, 400)
Y = np.linspace(0, 600, 400)
# X-nullcline: dX/dt = 0  ->  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX
Xnull = (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX
# Y-nullcline: dY/dt = 0  ->  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY
Ynull = (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

plt.figure(figsize=(7, 6))
plt.plot(Xnull, Y, 'b-', label='dX/dt = 0 nullcline')
plt.plot(X, Ynull, 'r-', label='dY/dt = 0 nullcline')
for s, eig, stable in results:
    color = 'green' if stable else 'black'
    marker = 'o' if stable else 'x'
    plt.plot(s[0], s[1], marker=marker, color=color, markersize=12,
             markeredgewidth=2,
             label=('stable' if stable else 'unstable (saddle)'))
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Toggle switch: nullclines and steady-state stability')
plt.legend()
plt.xlim(0, 600)
plt.ylim(0, 600)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3B.3.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print()
print("Explanation: The check confirms the result because a saddle point "
      "(one positive eigenvalue) repels trajectories along its unstable "
      "direction toward the two stable nodes, which is exactly the "
      "bistable behavior a toggle switch exhibits in simulation.")
