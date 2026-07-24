import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: self-activating gene ----
# basal transcription (g0) + excitatory Hill activation - linear degradation (k*X)
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4.0):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X

# ---- False-position (regula falsi) root finder on a single bracket ----
# Requires a sign change on [xmin, xmax]. Like bisection, but instead of the
# midpoint we use the x where the secant line through the endpoints hits zero.
def false_position(func, xmin, xmax, tol=1e-8, max_iter=200):
    f1 = func(xmin)          # value at left endpoint
    f2 = func(xmax)          # value at right endpoint
    if f1 * f2 > 0:          # no guaranteed root: endpoints share sign
        return None, 0
    it = 0
    x_new = xmin
    while it < max_iter:
        it += 1
        # secant/false-position update: zero crossing of line through endpoints
        x_new = (xmin * f2 - xmax * f1) / (f2 - f1)
        f_new = func(x_new)
        if abs(f_new) < tol:                # converged on the function value
            break
        if f1 * f_new < 0:                  # root is in [xmin, x_new]
            xmax, f2 = x_new, f_new
        else:                               # root is in [x_new, xmax]
            xmin, f1 = x_new, f_new
    return x_new, it

# ---- Windowed scan: split [a,b] into windows, run false position where a sign change exists ----
def windowed_scan(func, a, b, n_windows=200, tol=1e-8):
    edges = np.linspace(a, b, n_windows + 1)
    roots = []
    total_it = 0
    for i in range(n_windows):
        lo, hi = edges[i], edges[i + 1]
        if func(lo) * func(hi) <= 0:        # a sign change lives in this window
            r, it = false_position(func, lo, hi, tol=tol)
            total_it += it
            if r is not None:
                # avoid recording the same root twice at shared window edges
                if not any(abs(r - rr) < 1e-5 for rr in roots):
                    roots.append(r)
    return sorted(roots), total_it

k = 0.15
a, b = 0.0, 600.0
func = lambda X: f(X, k)

# ---- Whole-interval false position (single bracket [0, 600]) ----
whole_root, whole_it = false_position(func, a, b)
print(f"Whole-interval false position root: {whole_root}")
print(f"Whole-interval iterations: {whole_it}")
print(f"f(whole_root) = {func(whole_root):.3e}")

# ---- Windowed scan finds all roots ----
scan_roots, scan_it = windowed_scan(func, a, b)
print(f"Windowed-scan number of roots found: {len(scan_roots)}")
for j, r in enumerate(scan_roots):
    print(f"Windowed-scan root {j+1}: {r:.10f}   f = {func(r):.3e}")
print(f"Windowed-scan total iterations: {scan_it}")

# ---- Check: whole interval returns only ONE of the three roots ----
print(f"Check - roots on whole interval (single bracket): 1")
print(f"Check - roots via windowed scan: {len(scan_roots)}")
whole_matches_index = min(range(len(scan_roots)),
                          key=lambda j: abs(scan_roots[j] - whole_root))
print(f"Whole-interval root matches windowed root #{whole_matches_index+1}")
# One sentence explanation of why the check confirms the result:
print("Explanation: A single bracket over [0,600] has one net sign change, so "
      "false position can only converge to one root, while splitting the interval "
      "into windows isolates each individual sign change and recovers all three.")

# ---- Plot ----
X = np.linspace(a, b, 2000)
plt.figure(figsize=(9, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, func(X), label="f(X, k=0.15)")
plt.plot(scan_roots, [func(r) for r in scan_roots], "go", ms=9, label="windowed roots")
plt.plot([whole_root], [func(whole_root)], "rx", ms=12, mew=3, label="whole-interval root")
plt.xlabel("X")
plt.ylabel("f(X)")
plt.title("Self-activating gene: false-position roots at k=0.15")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.4.1_s2.png")
