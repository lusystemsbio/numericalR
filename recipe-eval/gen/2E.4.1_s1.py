import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: self-activating gene ---
# f(X, k) = g0 + g1*(X/Xth)^n/(1 + (X/Xth)^n) - k*X
g0, g1, Xth, n = 10.0, 45.0, 200.0, 4
k = 0.15

def f(X, k):
    r = (X / Xth) ** n
    return g0 + g1 * r / (1.0 + r) - k * X

# --- False-position (regula falsi) on a single bracket [xmin, xmax] ---
# Like bisection, but the test point is where the straight line through
# the two endpoints crosses zero, instead of the midpoint.
def false_position(func, xmin, xmax, tol=1e-8, max_iter=200):
    f1 = func(xmin)          # value at left end
    f2 = func(xmax)          # value at right end
    if f1 * f2 > 0:          # no sign change -> no guaranteed root here
        return None, 0
    x_new = xmin
    for it in range(1, max_iter + 1):
        # linear-interpolation crossing point (this is the whole method)
        x_new = (xmin * f2 - xmax * f1) / (f2 - f1)
        f_new = func(x_new)
        if abs(f_new) < tol:            # converged on the root value
            return x_new, it
        # keep the sub-interval that still brackets the root
        if f1 * f_new < 0:
            xmax, f2 = x_new, f_new     # root is in [xmin, x_new]
        else:
            xmin, f1 = x_new, f_new     # root is in [x_new, xmax]
    return x_new, max_iter

# --- Windowed scan: split [a,b] into windows, run false position where sign flips ---
def windowed_scan(func, a, b, num_windows=60, tol=1e-8):
    edges = np.linspace(a, b, num_windows + 1)
    roots, iters = [], []
    for i in range(num_windows):
        lo, hi = edges[i], edges[i + 1]
        if func(lo) * func(hi) <= 0:        # sign change -> a root lives here
            r, it = false_position(func, lo, hi, tol=tol)
            if r is not None:
                # avoid recording the same root twice at a shared window edge
                if not any(abs(r - rr) < 1e-6 for rr in roots):
                    roots.append(r)
                    iters.append(it)
    return roots, iters

# --- Run at k = 0.15 ---
# (1) whole interval, single bracket
root_whole, iters_whole = false_position(f, 0.0, 600.0)
print(f"Whole-interval false position root: {root_whole:.6f}  (iterations: {iters_whole})")
print(f"  f(root) = {f(root_whole, k):.3e}")

# (2) windowed scan
roots_win, iters_win = windowed_scan(f, 0.0, 600.0)
print(f"Number of roots found by windowed scan: {len(roots_win)}")
for j, (r, it) in enumerate(zip(roots_win, iters_win), 1):
    print(f"  windowed root {j}: {r:.6f}  (iterations: {it}, f = {f(r, k):.3e})")

# --- Check: whole-interval returns exactly ONE of the three roots ---
whole_is_subset = any(abs(root_whole - r) < 1e-4 for r in roots_win)
print(f"Whole-interval found {1} root; windowed scan found {len(roots_win)} roots.")
print(f"Check passes (single root is one of the windowed roots): {whole_is_subset and len(roots_win) == 3}")

# One-sentence explanation:
print("Explanation: The check confirms the result because a single sign-change "
      "bracketing method can only converge to one root in [0,600], while the "
      "windowed scan isolates each sign change separately and recovers all three, "
      "so agreement of the lone whole-interval root with one member of the "
      "three-root set proves the windowing is what exposes the extra roots.")

# --- Plot ---
X = np.linspace(0, 600, 1000)
plt.figure(figsize=(9, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X, k), label="f(X, k=0.15)")
plt.plot(root_whole, 0, "s", ms=11, mfc="none", mec="red",
         label=f"whole-interval root ({root_whole:.1f})")
plt.plot(roots_win, [0] * len(roots_win), "o", color="black",
         label="windowed roots")
plt.xlabel("X")
plt.ylabel("f(X, k)")
plt.title("False-position roots: whole interval vs windowed scan (k=0.15)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.4.1_s1.png")
