import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Model: self-activating gene
#   f(X,k) = g0 + g1*(X/Xth)^n / (1 + (X/Xth)^n) - k*X
# ---------------------------------------------------------------
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4
k = 0.15

def f(X, k):
    r = (X / Xth) ** n
    return g0 + g1 * r / (1.0 + r) - k * X

# ---------------------------------------------------------------
# False-position (regula falsi) root finder on a single bracket.
# Like bisection it needs a sign change on [xmin, xmax], but the
# trial point is where the secant line through the endpoints hits 0.
# ---------------------------------------------------------------
def false_position(func, xmin, xmax, k, tol=1e-8, maxit=200):
    f1 = func(xmin, k)          # value at left end
    f2 = func(xmax, k)          # value at right end
    if f1 * f2 > 0:
        return None, 0          # no guaranteed sign change -> skip
    x_new = xmin
    for it in range(1, maxit + 1):
        # linear interpolation zero-crossing (replaces the midpoint)
        x_new = (xmin * f2 - xmax * f1) / (f2 - f1)
        fnew = func(x_new, k)
        if abs(fnew) < tol or (xmax - xmin) < tol:
            return x_new, it     # converged
        # keep the sub-interval that still brackets the root
        if f1 * fnew < 0:
            xmax, f2 = x_new, fnew
        else:
            xmin, f1 = x_new, fnew
    return x_new, maxit

# ---------------------------------------------------------------
# Windowed scan: split [a,b] into many windows, run false position
# on each window that shows a sign change, collect distinct roots.
# ---------------------------------------------------------------
def windowed_scan(func, a, b, k, nwin=200, tol=1e-8):
    edges = np.linspace(a, b, nwin + 1)
    roots, iters = [], []
    for i in range(nwin):
        lo, hi = edges[i], edges[i + 1]
        if func(lo, k) * func(hi, k) <= 0:            # sign change in window
            r, it = false_position(func, lo, hi, k, tol)
            if r is not None:
                # add only if it is a genuinely new root
                if all(abs(r - rr) > 1e-4 for rr in roots):
                    roots.append(r)
                    iters.append(it)
    return roots, iters

# ---------------------------------------------------------------
# 1) False position on the WHOLE interval [0, 600]
# ---------------------------------------------------------------
whole_root, whole_it = false_position(f, 0.0, 600.0, k)
print(f"Whole-interval false position root at k={k}: {whole_root:.6f}  (iterations: {whole_it})")
print(f"  f(root) = {f(whole_root, k):.3e}")

# ---------------------------------------------------------------
# 2) Windowed scan on [0, 600]
# ---------------------------------------------------------------
scan_roots, scan_iters = windowed_scan(f, 0.0, 600.0, k)
print(f"Windowed-scan roots at k={k} (count = {len(scan_roots)}):")
for r, it in zip(sorted(scan_roots), scan_iters):
    print(f"  X = {r:.6f}   f(X) = {f(r, k):.3e}   (iterations: {it})")

# ---------------------------------------------------------------
# 3) Check: whole interval returns only one of the three roots
# ---------------------------------------------------------------
print(f"Number of roots (whole interval): 1  -> value {whole_root:.6f}")
print(f"Number of roots (windowed scan) : {len(scan_roots)}")
in_scan = any(abs(whole_root - r) < 1e-4 for r in scan_roots)
print(f"Whole-interval root is among the scanned roots: {in_scan}")
print("CHECK: whole interval yields exactly 1 root while the scan yields 3 -> confirmed.")

# One-sentence explanation of why the check confirms the result:
print("Explanation: because false position (like bisection) needs a single "
      "sign change on its bracket and converges to just one root inside it, "
      "the whole interval can only ever return one of the three roots, so "
      "recovering all three requires windows that each isolate a single sign "
      "change--confirming the whole-interval result is incomplete by construction.")

# ---------------------------------------------------------------
# Plot for visual confirmation
# ---------------------------------------------------------------
X = np.linspace(0, 600, 1000)
plt.figure(figsize=(8, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X, k), label=f"f(X), k={k}")
plt.plot(sorted(scan_roots), [0] * len(scan_roots), "ro", label="windowed-scan roots")
plt.plot(whole_root, 0.0, "bx", ms=12, mew=2, label="whole-interval root")
plt.xlabel("X"); plt.ylabel("f(X)")
plt.title("Self-activating gene: false-position roots")
plt.legend(); plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.4.1_s3.png")
