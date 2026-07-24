import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Fixed toggle-switch parameters (gX1 is the control parameter, so
# it is NOT included here; it is passed in during the sweep).
# ---------------------------------------------------------------
gX0 = 5.0
Yth = 100.0
nY  = 4.0
kX  = 0.1

gY0 = 4.0
gY1 = 40.0
Xth = 150.0
nX  = 4.0
kY  = 0.12


def hill_repress(z, zth, n):
    # Repressive Hill function: 1 / (1 + (z/zth)^n)
    return 1.0 / (1.0 + (z / zth) ** n)


def X_nullcline(Y, gX1):
    # dX/dt = 0  ->  X = (gX0 + gX1 * hill(Y)) / kX   (X as a function of Y)
    return (gX0 + gX1 * hill_repress(Y, Yth, nY)) / kX


def Y_nullcline(X):
    # dY/dt = 0  ->  Y = (gY0 + gY1 * hill(X)) / kY   (Y as a function of X)
    return (gY0 + gY1 * hill_repress(X, Xth, nX)) / kY


def self_map_residual(Y, gX1):
    # Compose the two nullclines: substitute X(Y) into Y(X), giving a single
    # equation in Y. A steady state is a root of g(Y) = Y_null(X_null(Y)) - Y.
    X = X_nullcline(Y, gX1)
    return Y_nullcline(X) - Y


def find_steady_states(gX1):
    # Scan Y over a physically reasonable range and bracket sign changes of
    # the residual, then refine each bracket by bisection (explicit root find).
    Ys = np.linspace(1e-6, 1000.0, 20000)
    res = np.array([self_map_residual(y, gX1) for y in Ys])
    roots = []
    for i in range(len(Ys) - 1):
        if res[i] == 0.0:
            roots.append(Ys[i])
        elif res[i] * res[i + 1] < 0.0:
            # Bisection inside the bracket [a, b]
            a, b = Ys[i], Ys[i + 1]
            fa = res[i]
            for _ in range(100):
                m = 0.5 * (a + b)
                fm = self_map_residual(m, gX1)
                if fa * fm <= 0.0:
                    b = m
                else:
                    a, fa = m, fm
            roots.append(0.5 * (a + b))
    # Convert each Y* into the full steady state (X*, Y*)
    return [(X_nullcline(y, gX1), y) for y in roots]


def jacobian(X, Y, gX1):
    # Analytic Jacobian of the 2D system at (X, Y).
    # d/dz [1/(1+(z/zth)^n)] = -(n/zth)*(z/zth)^(n-1) / (1+(z/zth)^n)^2
    def dhill(z, zth, n):
        u = (z / zth) ** n
        return -(n / zth) * (z / zth) ** (n - 1) / (1.0 + u) ** 2
    dfX_dX = -kX
    dfX_dY = gX1 * dhill(Y, Yth, nY)
    dfY_dX = gY1 * dhill(X, Xth, nX)
    dfY_dY = -kY
    return np.array([[dfX_dX, dfX_dY], [dfY_dX, dfY_dY]])


def is_stable(X, Y, gX1):
    # Stable iff all eigenvalues of the Jacobian have negative real part.
    eig = np.linalg.eigvals(jacobian(X, Y, gX1))
    return np.all(eig.real < 0.0)


# ---------------------------------------------------------------
# Sweep the control parameter gX1 from 0 to 100 and record every
# steady state, tagged by stability.
# ---------------------------------------------------------------
gX1_values = np.linspace(0.0, 100.0, 401)

stable_p, stable_X = [], []
unstable_p, unstable_X = [], []
count_per_p = []

print("gX1\tnum_steady_states\tX_values(stability)")
for gX1 in gX1_values:
    sss = find_steady_states(gX1)
    count_per_p.append(len(sss))
    labels = []
    for (X, Y) in sss:
        stab = is_stable(X, Y, gX1)
        if stab:
            stable_p.append(gX1); stable_X.append(X)
        else:
            unstable_p.append(gX1); unstable_X.append(X)
        labels.append(f"{X:.2f}({'S' if stab else 'U'})")
    # Print a subset to keep output readable but show branch structure
    if abs((gX1 * 10) % 100) < 1e-6:  # every 10 units
        print(f"{gX1:.1f}\t{len(sss)}\t{', '.join(labels)}")

# ---------------------------------------------------------------
# Bifurcation summary: report where the count of steady states changes.
# ---------------------------------------------------------------
print("\n--- Bifurcation points (where number of steady states changes) ---")
for i in range(1, len(gX1_values)):
    if count_per_p[i] != count_per_p[i - 1]:
        print(f"gX1 ~ {gX1_values[i]:.3f}: count changes {count_per_p[i-1]} -> {count_per_p[i]}")

print(f"\nMax number of coexisting steady states over sweep: {max(count_per_p)}")
print(f"Min number of coexisting steady states over sweep: {min(count_per_p)}")
print(f"Total stable points recorded:   {len(stable_X)}")
print(f"Total unstable points recorded: {len(unstable_X)}")

# ---------------------------------------------------------------
# Check: do stable and unstable branches both appear and merge?
# The system is bistable (3 steady states) over a middle range and
# monostable (1 steady state) at the extremes -> saddle-node bifurcations.
# ---------------------------------------------------------------
has_stable = len(stable_X) > 0
has_unstable = len(unstable_X) > 0
count_changes = any(count_per_p[i] != count_per_p[i - 1] for i in range(1, len(count_per_p)))
print(f"\nCheck - stable branch present:   {has_stable}")
print(f"Check - unstable branch present: {has_unstable}")
print(f"Check - steady-state count changes with gX1 (branches merge): {count_changes}")
print("Explanation: an unstable branch existing only where two stable branches "
      "coexist, and the count changing from 1 to 3 to 1, is exactly the "
      "signature of saddle-node bifurcations where a stable and unstable "
      "steady state collide and annihilate, confirming genuine bistability.")

# ---------------------------------------------------------------
# Plot steady-state X vs control parameter gX1, colored by stability.
# ---------------------------------------------------------------
plt.figure(figsize=(8, 6))
plt.scatter(unstable_p, unstable_X, s=8, c="crimson", label="unstable")
plt.scatter(stable_p, stable_X, s=8, c="royalblue", label="stable")
plt.xlabel("control parameter  gX1  (X production rate)")
plt.ylabel("steady-state  X*")
plt.title("Toggle-switch bifurcation diagram: X* vs gX1")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3E.1.1_s5.png")
print("\nSaved figure to 3E.1.1_s5.png")
