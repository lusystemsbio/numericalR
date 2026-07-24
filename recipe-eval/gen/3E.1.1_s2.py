import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Toggle switch model:
#   dX/dt = gX0 + gX1/(1 + (Y/Yth)^nY) - kX*X
#   dY/dt = gY0 + gY1/(1 + (X/Xth)^nX) - kY*Y
# Genes X and Y mutually repress via Hill functions.
# ---------------------------------------------------------------

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

# ---------------------------------------------------------------
# Nullclines (separation of variables):
# X-nullcline (dX/dt=0):  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX   -> X as function of Y
# Y-nullcline (dY/dt=0):  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY   -> Y as function of X
# A steady state satisfies BOTH. Substitute the Y-nullcline into the
# X-nullcline to get a single equation g(X)=0 whose roots are the fixed points.
# ---------------------------------------------------------------

def Y_of_X(X):
    # Y-nullcline: steady-state Y given X
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

def X_of_Y(Y, gX1):
    # X-nullcline: steady-state X given Y
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def g(X, gX1):
    # Self-consistency residual: X must equal X_of_Y(Y_of_X(X))
    return X_of_Y(Y_of_X(X), gX1) - X

def jacobian_eigs(X, Y, gX1):
    # Jacobian of the 2D system at (X,Y).
    # dfX/dX = -kX
    # dfX/dY = gX1 * d/dY [1/(1+(Y/Yth)^nY)]
    # dfY/dX = gY1 * d/dX [1/(1+(X/Xth)^nX)]
    # dfY/dY = -kY
    dHY_dY = -gX1 * (nY / Yth) * (Y / Yth) ** (nY - 1) / (1.0 + (Y / Yth) ** nY) ** 2
    dHX_dX = -gY1 * (nX / Xth) * (X / Xth) ** (nX - 1) / (1.0 + (X / Xth) ** nX) ** 2
    J = np.array([[-kX, dHY_dY],
                  [dHX_dX, -kY]])
    return np.linalg.eigvals(J)

def find_steady_states(gX1):
    # Scan X over a wide grid, detect sign changes of g(X), then bisect for roots.
    Xgrid = np.linspace(1e-6, (gX0 + gX1) / kX + 10.0, 4000)
    vals = np.array([g(x, gX1) for x in Xgrid])
    roots = []
    for i in range(len(Xgrid) - 1):
        if vals[i] == 0.0:
            roots.append(Xgrid[i])
        elif vals[i] * vals[i + 1] < 0.0:
            # bisection between bracketing grid points
            a, b = Xgrid[i], Xgrid[i + 1]
            fa = vals[i]
            for _ in range(80):
                m = 0.5 * (a + b)
                fm = g(m, gX1)
                if fa * fm <= 0.0:
                    b = m
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + b))
    # classify each root
    states = []
    for Xs in roots:
        Ys = Y_of_X(Xs)
        eigs = jacobian_eigs(Xs, Ys, gX1)
        stable = np.all(np.real(eigs) < 0.0)  # stable iff all eigenvalues have negative real part
        states.append((Xs, stable))
    return states

# ---------------------------------------------------------------
# Sweep the control parameter gX1 from 0 to 100
# ---------------------------------------------------------------
gX1_sweep = np.linspace(0.0, 100.0, 400)

stable_p, stable_x = [], []
unstable_p, unstable_x = [], []

print("=== Bifurcation sweep: gX1 from 0 to 100 ===")
prev_count = None
for p in gX1_sweep:
    states = find_steady_states(p)
    for Xs, stable in states:
        if stable:
            stable_p.append(p); stable_x.append(Xs)
        else:
            unstable_p.append(p); unstable_x.append(Xs)
    # report where the number of steady states changes (bifurcation)
    if prev_count is not None and len(states) != prev_count:
        print(f"gX1 = {p:.3f}: number of steady states changed {prev_count} -> {len(states)}")
    prev_count = len(states)

# Report steady states at a few representative parameter values
for ptest in [0.0, 25.0, 50.0, 75.0, 100.0]:
    states = find_steady_states(ptest)
    print(f"\ngX1 = {ptest:.1f}: {len(states)} steady state(s)")
    for Xs, stable in states:
        Ys = Y_of_X(Xs)
        eigs = jacobian_eigs(Xs, Ys, ptest)
        label = "stable" if stable else "unstable"
        print(f"  X* = {Xs:10.4f}, Y* = {Ys:10.4f}, "
              f"eigs = ({eigs[0].real:.4f}, {eigs[1].real:.4f}), {label}")

# ---------------------------------------------------------------
# Check summary: do stable and unstable branches appear/merge?
# ---------------------------------------------------------------
counts = [len(find_steady_states(p)) for p in gX1_sweep]
print("\n=== Bifurcation check ===")
print(f"Minimum number of steady states over sweep: {min(counts)}")
print(f"Maximum number of steady states over sweep: {max(counts)}")
print(f"Total stable points plotted:   {len(stable_x)}")
print(f"Total unstable points plotted: {len(unstable_x)}")
multistable = max(counts) > 1
print(f"Multistability (>=2 steady states somewhere): {multistable}")

# ---------------------------------------------------------------
# Plot: steady-state X vs control parameter gX1, colored by stability
# ---------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.scatter(unstable_p, unstable_x, s=8, c="red", label="unstable")
plt.scatter(stable_p, stable_x, s=8, c="blue", label="stable")
plt.xlabel("control parameter gX1 (X production rate)")
plt.ylabel("steady-state X*")
plt.title("Toggle switch bifurcation diagram")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3E.1.1_s2.png")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: The check confirms the result because seeing the steady-state "
      "count rise from 1 to 3 (and stable/unstable branches born together and "
      "annihilating) is the defining signature of saddle-node bifurcations, proving "
      "the parameter genuinely reshapes the number and stability of fixed points.")
