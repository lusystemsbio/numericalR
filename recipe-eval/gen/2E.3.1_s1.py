import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: self-activating gene ----
# f(X, k) = basal + excitatory Hill activation - linear degradation
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4.0):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X

# ---- Homemade bisection root-finder ----
# Given an interval [a, b] where f changes sign, repeatedly halve it,
# keeping the half that still brackets the root, until it is tiny.
def bisection(func, a, b, k, tol=1e-10, maxit=200):
    fa = func(a, k)
    fb = func(b, k)
    if fa == 0.0:  # endpoint is exactly a root
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:  # no sign change -> cannot bracket a root here
        return None
    for _ in range(maxit):
        m = 0.5 * (a + b)      # midpoint
        fm = func(m, k)
        if fm == 0.0 or (b - a) < tol:  # converged
            return m
        if fa * fm < 0.0:      # root is in the left half
            b, fb = m, fm
        else:                  # root is in the right half
            a, fa = m, fm
    return 0.5 * (a + b)

# ---- Whole-interval bisection: treat [lo, hi] as a single bracket ----
def whole_interval_root(func, lo, hi, k):
    return bisection(func, lo, hi, k)  # returns at most one root

# ---- Windowed scan: split [lo, hi] into small windows and bisect ----
# each window that shows a sign change, recovering every root.
def windowed_roots(func, lo, hi, k, nwin=600, tol=1e-10):
    edges = np.linspace(lo, hi, nwin + 1)  # window boundaries
    roots = []
    for i in range(nwin):
        a, b = edges[i], edges[i + 1]
        if func(a, k) * func(b, k) <= 0.0:  # sign change in this window
            r = bisection(func, a, b, k, tol=tol)
            if r is not None:
                # avoid duplicates at shared window edges
                if not any(abs(r - rr) < 1e-6 for rr in roots):
                    roots.append(r)
    return sorted(roots)

lo, hi = 0.0, 600.0
ks = [0.12, 0.15, 0.2]

results = {}
for k in ks:
    whole = whole_interval_root(f, lo, hi, k)
    scan = windowed_roots(f, lo, hi, k)
    results[k] = (whole, scan)
    print(f"k = {k}")
    print(f"  whole-interval bisection root(s): {[round(whole, 6)] if whole is not None else []}")
    print(f"  windowed-scan roots ({len(scan)} found): {[round(r, 6) for r in scan]}")

# ---- Explicit check ----
print()
print("CHECK:")
for k in [0.12, 0.2]:
    whole, scan = results[k]
    print(f"  monostable k = {k}: whole-interval finds 1 root ({round(whole,6)}), "
          f"scan finds {len(scan)} root(s) -> agree: {len(scan) == 1}")
k = 0.15
whole, scan = results[k]
print(f"  bistable k = {k}: whole-interval finds only 1 of the roots ({round(whole,6)}), "
      f"scan finds {len(scan)} roots: {[round(r,6) for r in scan]}")

# Explanation (one sentence):
print()
print("Explanation: whole-interval bisection can only ever return one root because a single")
print("sign change over [0,600] hides the even number of interior crossings at the bistable k,")
print("so the fact that the windowed scan recovers three roots at k=0.15 (vs one at k=0.12,0.2)")
print("confirms that k=0.15 is genuinely bistable and that the scan captures roots bisection alone misses.")

# ---- Plot ----
X = np.linspace(lo, hi, 1000)
fig, ax = plt.subplots(figsize=(9, 6))
colors = {0.12: "tab:blue", 0.15: "tab:red", 0.2: "tab:green"}
for k in ks:
    ax.plot(X, f(X, k), color=colors[k], label=f"f(X), k={k}")
    whole, scan = results[k]
    if whole is not None:
        ax.plot(whole, 0.0, "o", color=colors[k], mfc="none", ms=13,
                label=f"whole-interval root (k={k})")
    ax.plot(scan, np.zeros(len(scan)), "x", color=colors[k], ms=10,
            label=f"scan roots (k={k})")
ax.axhline(0.0, color="k", lw=0.8)
ax.set_xlabel("X")
ax.set_ylabel("f(X, k)")
ax.set_title("Self-activating gene steady states: whole-interval vs windowed bisection")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.3.1_s1.png")
