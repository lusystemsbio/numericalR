import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# --- Model: self-activating gene ---
# f(X,k) = basal g0 + excitatory Hill term - linear degradation k*X
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4):
    r = (X / Xth) ** n
    return g0 + g1 * r / (1.0 + r) - k * X


# --- Homemade bisection root-finder ---
# Requires f(a) and f(b) to have opposite signs (a sign change brackets a root).
def bisect(func, a, b, k, tol=1e-9, maxit=200):
    fa = func(a, k)
    fb = func(b, k)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:
        return None  # no sign change -> not bracketed, bail out
    for _ in range(maxit):
        m = 0.5 * (a + b)          # midpoint
        fm = func(m, k)
        if fm == 0.0 or 0.5 * (b - a) < tol:
            return m               # converged
        # keep the half that still brackets the root (sign change with fa)
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


# --- Whole-interval solve: one bisection call on the full bracket if it changes sign ---
def whole_interval(func, a, b, k):
    root = bisect(func, a, b, k)
    return [] if root is None else [root]


# --- Windowed scan: split [a,b] into small windows, bisect each that changes sign ---
def windowed_scan(func, a, b, k, nwin=600):
    edges = np.linspace(a, b, nwin + 1)
    roots = []
    for i in range(nwin):
        lo, hi = edges[i], edges[i + 1]
        if func(lo, k) * func(hi, k) < 0.0:   # sign change inside this window
            r = bisect(func, lo, hi, k)
            if r is not None:
                # avoid near-duplicate roots at shared window edges
                if not any(abs(r - q) < 1e-6 for q in roots):
                    roots.append(r)
    return sorted(roots)


a, b = 0.0, 600.0
ks = [0.12, 0.15, 0.2]

results = {}
for k in ks:
    whole = whole_interval(f, a, b, k)
    scan = windowed_scan(f, a, b, k)
    results[k] = (whole, scan)
    print(f"k = {k}")
    print(f"  whole-interval bisection roots ({len(whole)}): "
          + ", ".join(f"{r:.6f}" for r in whole))
    print(f"  windowed-scan roots ({len(scan)}): "
          + ", ".join(f"{r:.6f}" for r in scan))

# --- Separate check ---
print()
print("CHECK:")
for k in ks:
    whole, scan = results[k]
    label = "bistable" if len(scan) == 3 else "monostable"
    print(f"  k = {k} ({label}): whole-interval found {len(whole)} root(s), "
          f"windowed scan found {len(scan)} root(s)")

print()
print("Explanation: The check confirms the result because at monostable k the "
      "whole-interval and windowed counts agree (1 root each), while at bistable "
      "k=0.15 the single whole-interval solve returns only one of the three roots "
      "(bisection converges to just one bracketed crossing) yet the windowed scan "
      "isolates each sign change separately and recovers all three.")

# --- Plot ---
fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
Xs = np.linspace(a, b, 1000)
for ax, k in zip(axes, ks):
    ax.axhline(0.0, color="gray", lw=0.8)
    ax.plot(Xs, [f(X, k) for X in Xs], label=f"f(X, k={k})")
    whole, scan = results[k]
    ax.scatter(scan, [0.0] * len(scan), s=120, facecolors="none",
               edgecolors="green", linewidths=2, label="windowed scan", zorder=5)
    ax.scatter(whole, [0.0] * len(whole), s=40, color="red",
               label="whole interval", zorder=6)
    ax.set_title(f"k = {k}  ({'bistable' if len(scan)==3 else 'monostable'})")
    ax.set_xlabel("X")
    ax.legend(fontsize=8)
axes[0].set_ylabel("f(X, k)")
fig.suptitle("Self-activating gene: steady states via homemade bisection")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.3.1_s4.png")
