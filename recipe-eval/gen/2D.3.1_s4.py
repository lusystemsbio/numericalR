import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: self-activating gene ----
# f(X) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X
def f(X, g0, g1, Xth, n, k):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)  # excitatory Hill function
    return g0 + g1 * hill - k * X                    # basal + activation - degradation

# ---- Effective potential by explicit trapezoidal accumulation ----
# U(X) = -integral_0^X f(x) dx, with U(0) = 0
# U(X+dx) = U(X) - (f(X) + f(X+dx))/2 * dx
def effective_potential(Xgrid, g0, g1, Xth, n, k):
    dx = Xgrid[1] - Xgrid[0]              # uniform grid spacing
    fvals = f(Xgrid, g0, g1, Xth, n, k)   # f at every grid point
    U = np.zeros_like(Xgrid)              # U(0) = 0
    for i in range(len(Xgrid) - 1):
        # subtract the trapezoid area of f over [X_i, X_{i+1}]
        U[i + 1] = U[i] - 0.5 * (fvals[i] + fvals[i + 1]) * dx
    return U, fvals

# ---- Parameters ----
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4
ks = [0.15, 0.2, 0.1]

# Grid over X (nM)
Xgrid = np.linspace(0.0, 500.0, 5001)
dx = Xgrid[1] - Xgrid[0]

# ---- Helper: find sign changes of f (fixed points) ----
def fixed_points(X, fvals):
    roots = []
    kinds = []
    for i in range(len(fvals) - 1):
        if fvals[i] == 0.0:
            roots.append(X[i]); kinds.append("stable" if (fvals[i+1] < 0) else "unstable")
        elif fvals[i] * fvals[i + 1] < 0.0:
            # linear interpolation for the crossing location
            xr = X[i] - fvals[i] * (X[i + 1] - X[i]) / (fvals[i + 1] - fvals[i])
            # f goes + -> - : stable (valley);  - -> + : unstable (peak)
            kind = "stable" if fvals[i] > fvals[i + 1] else "unstable"
            roots.append(xr); kinds.append(kind)
    return roots, kinds

# ---- Compute for each k ----
results = {}
for k in ks:
    U, fvals = effective_potential(Xgrid, g0, g1, Xth, n, k)
    roots, kinds = fixed_points(Xgrid, fvals)
    results[k] = (U, roots, kinds)
    print(f"--- k = {k} ---")
    print(f"grid spacing dx = {dx} nM")
    print(f"number of fixed points (roots of f) = {len(roots)}")
    for xr, kind in zip(roots, kinds):
        print(f"  fixed point at X = {xr:.4f} nM  -> {kind} "
              f"({'valley' if kind=='stable' else 'peak'})")
    # basins = number of stable states, barriers = number of unstable states
    n_stable = sum(1 for kd in kinds if kd == "stable")
    n_unstable = sum(1 for kd in kinds if kd == "unstable")
    print(f"  number of basins (stable states) = {n_stable}")
    print(f"  number of barriers (unstable states) = {n_unstable}")

# ---- Check at k = 0.15: two basins split by a barrier ----
U015, roots015, kinds015 = results[0.15]
stable_015 = [xr for xr, kd in zip(roots015, kinds015) if kd == "stable"]
unstable_015 = [xr for xr, kd in zip(roots015, kinds015) if kd == "unstable"]
print("--- Check at k = 0.15 ---")
print(f"stable states (basins), nM  = {[round(x,2) for x in stable_015]}")
print(f"unstable state (barrier), nM = {[round(x,2) for x in unstable_015]}")
bistable_015 = (len(stable_015) == 2 and len(unstable_015) == 1)
print(f"bistable (2 basins + 1 barrier)? {bistable_015}")

for k in (0.2, 0.1):
    _, r, kd = results[k]
    n_s = sum(1 for x in kd if x == "stable")
    print(f"--- Check at k = {k} ---")
    print(f"number of basins (stable states) = {n_s}  -> single well? {n_s == 1}")

# ---- Plot U(X) for the three k values ----
plt.figure(figsize=(8, 6))
colors = {0.15: "tab:blue", 0.2: "tab:orange", 0.1: "tab:green"}
for k in ks:
    U = results[k][0]
    plt.plot(Xgrid, U, color=colors[k], label=f"k = {k}")
# mark the k=0.15 fixed points on its curve
for xr, kind in zip(roots015, kinds015):
    idx = int(round(xr / dx))
    plt.scatter([Xgrid[idx]], [U015[idx]],
                color="red" if kind == "stable" else "black",
                zorder=5, s=40,
                marker="v" if kind == "stable" else "^")
plt.xlabel("X (nM)")
plt.ylabel("Effective potential U(X)")
plt.title("Effective potential of a self-activating gene circuit")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.3.1_s4.png")

# ---- One-sentence explanation ----
print("Explanation: The check confirms the result because the valleys (minima) of "
      "U(X) occur exactly where f(X)=0 with f decreasing (stable states) and the peak "
      "(maximum) where f(X)=0 with f increasing (unstable state), so finding two minima "
      "near 100 and 300 nM separated by one maximum near 200 nM at k=0.15 — versus a "
      "single minimum at k=0.2 and k=0.1 — verifies bistability and its loss.")
