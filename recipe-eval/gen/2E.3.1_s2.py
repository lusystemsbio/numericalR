import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Model: self-activating gene. Steady states satisfy f(X, k) = 0.
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4.0):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)  # excitatory Hill term
    return g0 + g1 * hill - k * X                    # + basal - linear degradation


# Homemade bisection: assumes f(a,k) and f(b,k) have opposite signs.
def bisect(a, b, k, tol=1e-8, maxit=200):
    fa = f(a, k)
    fb = f(b, k)
    if fa == 0.0:  # endpoint is exactly a root
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0.0:  # no sign change -> cannot guarantee a root
        return None
    for _ in range(maxit):
        m = 0.5 * (a + b)      # midpoint
        fm = f(m, k)
        if fm == 0.0 or 0.5 * (b - a) < tol:
            return m
        # keep the half that still brackets the root (sign change preserved)
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    return 0.5 * (a + b)


# Windowed scan: split [lo, hi] into many small windows, bisect each one
# that shows a sign change, so every root is recovered.
def windowed_roots(lo, hi, k, nwin=600, tol=1e-8):
    edges = np.linspace(lo, hi, nwin + 1)
    roots = []
    for i in range(nwin):
        a, b = edges[i], edges[i + 1]
        if f(a, k) * f(b, k) <= 0.0:  # sign change (or endpoint zero) in window
            r = bisect(a, b, k, tol)
            if r is not None:
                # avoid duplicates from shared window edges
                if not any(abs(r - rr) < 1e-5 for rr in roots):
                    roots.append(r)
    return sorted(roots)


lo, hi = 0.0, 600.0
ks = [0.12, 0.15, 0.2]

print("=== Roots on whole interval [0, 600] (single bisection) vs windowed scan ===")
results = {}
for k in ks:
    whole = bisect(lo, hi, k)              # one bisection on the entire interval
    windowed = windowed_roots(lo, hi, k)   # scan of small windows
    results[k] = (whole, windowed)
    print(f"k = {k}: whole-interval bisection root = "
          f"{('%.6f' % whole) if whole is not None else 'None'}")
    print(f"k = {k}: windowed-scan roots ({len(windowed)}) = "
          f"[{', '.join('%.6f' % r for r in windowed)}]")

print()
print("=== Check: monostable vs bistable ===")
for k in ks:
    whole, windowed = results[k]
    n_windowed = len(windowed)
    kind = "monostable" if n_windowed == 1 else f"bistable ({n_windowed} states)"
    print(f"k = {k}: {kind}; whole-interval bisection found "
          f"{0 if whole is None else 1} root, windowed found {n_windowed}")

print()
print("Explanation: A single bisection on the whole interval can only return one "
      "root because it needs a single sign change across the endpoints, so at "
      "bistable k=0.15 (three roots) it finds just one while the windowed scan "
      "isolates each sign change separately and recovers all three; at monostable "
      "k=0.12 and k=0.2 both methods agree on the one existing root, which confirms "
      "that the extra roots at k=0.15 are real and not artifacts.")

# Plot f(X,k) with the windowed roots marked.
X = np.linspace(lo, hi, 2000)
plt.figure(figsize=(9, 6))
colors = ['tab:blue', 'tab:orange', 'tab:green']
for k, c in zip(ks, colors):
    plt.plot(X, f(X, k), color=c, label=f"f(X), k={k}")
    for r in results[k][1]:
        plt.plot(r, 0.0, 'o', color=c, ms=9, mec='k')
plt.axhline(0.0, color='k', lw=0.8)
plt.xlabel("X")
plt.ylabel("f(X, k)")
plt.title("Self-activating gene: steady states via windowed bisection")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.3.1_s2.png")
print()
print("Figure saved to 2E.3.1_s2.png")
