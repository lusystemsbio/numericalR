import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# --- Model: self-activating gene ---
# f(X) = basal + excitatory Hill activation - linear degradation
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X


# --- Homemade bisection root-finder ---
# Requires f(a) and f(b) to have opposite signs (a sign change brackets a root).
# Repeatedly halve the interval, keeping whichever half still brackets the root.
def bisect(func, a, b, tol=1e-8, maxit=200):
    fa, fb = func(a), func(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return None  # no sign change -> not bracketed, refuse
    for _ in range(maxit):
        m = 0.5 * (a + b)          # midpoint
        fm = func(m)
        if fm == 0.0 or 0.5 * (b - a) < tol:
            return m               # converged
        if fa * fm < 0.0:          # root is in the left half [a, m]
            b, fb = m, fm
        else:                      # root is in the right half [m, b]
            a, fa = m, fm
    return 0.5 * (a + b)


# --- Windowed scan: split [lo, hi] into small windows, ---
# --- run bisection on every window that shows a sign change. ---
def windowed_scan(func, lo, hi, nwin=600, tol=1e-8):
    edges = np.linspace(lo, hi, nwin + 1)
    roots = []
    for i in range(nwin):
        a, b = edges[i], edges[i + 1]
        if func(a) * func(b) < 0.0:      # sign change in this window
            r = bisect(func, a, b, tol=tol)
            if r is not None:
                # de-duplicate roots shared across adjacent windows
                if not any(abs(r - rr) < 1e-6 for rr in roots):
                    roots.append(r)
    return sorted(roots)


lo, hi = 0.0, 600.0
ks = [0.12, 0.15, 0.2]

whole_roots = {}
scan_roots = {}
for k in ks:
    func = lambda X, k=k: f(X, k)
    # Bisection on the WHOLE interval: only sees the outer sign change,
    # so it can return at most one root even if several exist.
    whole = bisect(func, lo, hi)
    whole_roots[k] = whole
    # Windowed scan recovers every root by localizing each sign change.
    scan_roots[k] = windowed_scan(func, lo, hi)

# --- Print results ---
for k in ks:
    print(f"k = {k}: whole-interval bisection root = {whole_roots[k]}")
for k in ks:
    rs = ", ".join(f"{r:.6f}" for r in scan_roots[k])
    print(f"k = {k}: windowed-scan roots ({len(scan_roots[k])}) = [{rs}]")

# --- Separate check ---
for k in ks:
    n_scan = len(scan_roots[k])
    kind = "bistable (3 states)" if n_scan == 3 else "monostable (1 state)"
    print(f"CHECK k = {k}: {kind}; whole-interval found "
          f"{0 if whole_roots[k] is None else 1} root, scan found {n_scan}")

print("EXPLANATION: The check confirms the result because a single bracketing "
      "over [0,600] detects only the one outermost sign change (thus one root, "
      "even when three exist at bistable k=0.15), whereas subdividing into small "
      "windows isolates each individual sign change and so recovers all roots.")

# --- Plot ---
X = np.linspace(lo, hi, 2000)
fig, ax = plt.subplots(figsize=(9, 6))
colors = {0.12: "tab:blue", 0.15: "tab:red", 0.2: "tab:green"}
for k in ks:
    ax.plot(X, f(X, k), color=colors[k], label=f"f(X), k={k}")
    for r in scan_roots[k]:
        ax.plot(r, 0.0, "o", color=colors[k], ms=9,
                markeredgecolor="k", zorder=5)
    if whole_roots[k] is not None:
        ax.plot(whole_roots[k], 0.0, "x", color="black", ms=12, mew=2, zorder=6)
ax.axhline(0.0, color="gray", lw=0.8)
ax.set_xlabel("X")
ax.set_ylabel("f(X, k)")
ax.set_title("Self-activating gene: roots by whole-interval (x) vs windowed scan (o)")
ax.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.3.1_s5.png")
