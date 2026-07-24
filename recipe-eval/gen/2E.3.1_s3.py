import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: self-activating gene ---
# f(X, k) = basal + excitatory Hill activation - linear degradation
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4.0):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X

# --- Homemade bisection root-finder ---
# Requires a sign change on [a, b]; repeatedly halve, keeping the bracketing half.
def bisect(func, a, b, k, tol=1e-10, maxit=200):
    fa, fb = func(a, k), func(b, k)
    if fa == 0.0:            # endpoint is exactly a root
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:        # no sign change -> cannot bracket a root here
        return None
    for _ in range(maxit):
        m = 0.5 * (a + b)    # midpoint
        fm = func(m, k)
        if fm == 0.0 or 0.5 * (b - a) < tol:  # converged
            return m
        if fa * fm < 0.0:    # root lies in [a, m] -> keep left half
            b, fb = m, fm
        else:                # root lies in [m, b] -> keep right half
            a, fa = m, fm
    return 0.5 * (a + b)

# --- Whole-interval bisection: single attempt on [a, b] ---
def whole_interval(func, a, b, k):
    r = bisect(func, a, b, k)
    return [] if r is None else [r]

# --- Windowed scan: split [a, b] into windows, bisect each with a sign change ---
def windowed_scan(func, a, b, k, nwin=600, tol=1e-10):
    edges = np.linspace(a, b, nwin + 1)  # window boundaries
    roots = []
    for i in range(nwin):
        lo, hi = edges[i], edges[i + 1]
        if func(lo, k) * func(hi, k) <= 0.0:  # sign change in this window
            r = bisect(func, lo, hi, k)
            if r is not None:
                # avoid duplicates at shared window boundaries
                if not any(abs(r - q) < 1e-6 for q in roots):
                    roots.append(r)
    return sorted(roots)

# --- Run the tests ---
a, b = 0.0, 600.0
for k in (0.12, 0.15, 0.20):
    whole = whole_interval(f, a, b, k)
    scan = windowed_scan(f, a, b, k)
    print(f"k = {k}")
    print(f"  whole-interval roots (count {len(whole)}): "
          + ", ".join(f"{r:.6f}" for r in whole))
    print(f"  windowed-scan  roots (count {len(scan)}): "
          + ", ".join(f"{r:.6f}" for r in scan))

# --- Separate check: monostable vs bistable ---
print("\nCheck (whole-interval count, windowed-scan count):")
for k in (0.12, 0.15, 0.20):
    whole = whole_interval(f, a, b, k)
    scan = windowed_scan(f, a, b, k)
    label = "bistable" if len(scan) == 3 else "monostable"
    print(f"  k = {k}: whole={len(whole)}, scan={len(scan)}  -> {label}")

# --- Plot ---
X = np.linspace(a, b, 1000)
fig, ax = plt.subplots(figsize=(8, 5))
for k, c in zip((0.12, 0.15, 0.20), ("tab:blue", "tab:orange", "tab:green")):
    ax.plot(X, f(X, k), color=c, label=f"f(X), k={k}")
    for r in windowed_scan(f, a, b, k):
        ax.plot(r, 0.0, "o", color=c, ms=8)
ax.axhline(0.0, color="k", lw=0.8)
ax.set_xlabel("X")
ax.set_ylabel("f(X, k)")
ax.set_title("Self-activating gene: steady states via bracketing bisection")
ax.legend()
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.3.1_s3.png")

# Explanation:
# The check confirms the result because a single bisection over [0,600] can only ever
# return one root (it commits to one bracketing half at each step and discards the
# other), so it silently misses the extra stable/unstable states at bistable k=0.15;
# the windowed scan isolates each sign change in its own sub-interval and thus recovers
# all three roots, proving the extra steady states are real and were merely hidden by
# the single-bracket method.
print("\nWhy the check confirms the result: a single whole-interval bisection can return "
      "only one of several roots because it keeps just one bracketing half at each step, "
      "so recovering all three roots at k=0.15 only via the windowed scan proves the "
      "extra steady states exist and were simply missed by the single-bracket search.")
