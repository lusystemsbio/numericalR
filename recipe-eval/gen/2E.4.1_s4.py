import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----- Model: self-activating gene -----
# f(X,k) = basal transcription + excitatory Hill activation - linear degradation
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4):
    hill = (X / Xth)**n / (1.0 + (X / Xth)**n)   # excitatory Hill function
    return g0 + g1 * hill - k * X                 # minus linear degradation k*X

# ----- False-position (regula falsi) bracketing method -----
# Like bisection, but instead of the midpoint we take the x where the straight
# line joining the two endpoints crosses zero.
def false_position(func, xmin, xmax, k, tol=1e-8, maxit=200):
    f1 = func(xmin, k)          # value at left endpoint
    f2 = func(xmax, k)          # value at right endpoint
    if f1 * f2 > 0:             # no sign change -> no bracketed root
        return None, 0
    x_new = xmin
    for it in range(1, maxit + 1):
        # linear-interpolation update: zero of the secant through the endpoints
        x_new = (xmin * f2 - xmax * f1) / (f2 - f1)
        f_new = func(x_new, k)
        if abs(f_new) < tol:    # converged on the function value
            return x_new, it
        # keep the sub-interval that still brackets the root (sign change)
        if f1 * f_new < 0:
            xmax, f2 = x_new, f_new
        else:
            xmin, f1 = x_new, f_new
    return x_new, maxit

# ----- Windowed scan: subdivide the interval, run false position on each bracket -----
def windowed_roots(func, a, b, k, nwin=60, tol=1e-8):
    edges = np.linspace(a, b, nwin + 1)   # window boundaries
    roots, iters = [], []
    for i in range(nwin):
        lo, hi = edges[i], edges[i + 1]
        if func(lo, k) * func(hi, k) <= 0:          # sign change inside window
            r, it = false_position(func, lo, hi, k, tol=tol)
            if r is not None:
                # avoid recording the same root twice from adjacent windows
                if not any(abs(r - rr) < 1e-4 for rr in roots):
                    roots.append(r)
                    iters.append(it)
    return roots, iters

k = 0.15
a, b = 0.0, 600.0

# (1) False position on the WHOLE interval
root_whole, it_whole = false_position(f, a, b, k)
print("=== False position on whole interval [0, 600], k = 0.15 ===")
if root_whole is not None:
    print(f"root (whole interval)        = {root_whole:.6f}   (iterations = {it_whole})")
else:
    print("no bracketed root on whole interval")

# (2) Windowed scan
roots_win, iters_win = windowed_roots(f, a, b, k)
print("\n=== Windowed scan on [0, 600], k = 0.15 ===")
for i, (r, it) in enumerate(zip(roots_win, iters_win)):
    print(f"root {i+1} (windowed scan)      = {r:.6f}   (iterations = {it})")
print(f"number of roots found (whole)    = {0 if root_whole is None else 1}")
print(f"number of roots found (windowed) = {len(roots_win)}")

# Confirmation check: is the single whole-interval root among the windowed roots?
in_set = any(abs(root_whole - r) < 1e-4 for r in roots_win)
print(f"\nwhole-interval root is one of the windowed roots: {in_set}")
# Explanation:
# The whole-interval call returns exactly ONE value that also appears in the
# windowed set of THREE, which shows regula falsi only finds a single root per
# bracket and needs the windowed scan to isolate every sign change separately.

# ----- Plot -----
X = np.linspace(a, b, 1000)
plt.figure(figsize=(8, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X, k), label="f(X, k=0.15)")
plt.plot(roots_win, [0]*len(roots_win), "ro", ms=9, label="windowed roots")
if root_whole is not None:
    plt.plot([root_whole], [0], "b*", ms=16, label="whole-interval root")
plt.xlabel("X"); plt.ylabel("f(X)")
plt.title("Self-activating gene: false-position roots (k = 0.15)")
plt.legend(); plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.4.1_s4.png")
