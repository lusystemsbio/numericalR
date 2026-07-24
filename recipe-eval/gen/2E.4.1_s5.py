import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: self-activating gene ---------------------------------------
# f(X,k) = basal transcription + excitatory Hill activation - linear degradation
def f(X, k, g0=10.0, g1=45.0, Xth=200.0, n=4.0):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)  # Hill activation term
    return g0 + g1 * hill - k * X                    # minus linear degradation k*X

# --- False-position (regula falsi) on a single bracket [xmin, xmax] ----
def false_position(func, xmin, xmax, k, tol=1e-8, maxit=200):
    f1 = func(xmin, k)                # value at left endpoint
    f2 = func(xmax, k)                # value at right endpoint
    if f1 * f2 > 0:                   # no guaranteed sign change -> no bracketed root
        return None, 0
    it = 0
    x_new = xmin
    for it in range(1, maxit + 1):
        # replace bisection midpoint with the secant-line zero crossing:
        x_new = (xmin * f2 - xmax * f1) / (f2 - f1)
        fn = func(x_new, k)
        if abs(fn) < tol:             # converged on function value
            break
        # keep the subinterval that still brackets the root (sign change)
        if f1 * fn < 0:
            xmax, f2 = x_new, fn      # root is in [xmin, x_new]
        else:
            xmin, f1 = x_new, fn      # root is in [x_new, xmax]
    return x_new, it

# --- Windowed driver: split domain into windows, bracket each sign change ---
def windowed_scan(func, a, b, k, nwin=200):
    roots, iters = [], []
    edges = np.linspace(a, b, nwin + 1)   # window boundaries
    for i in range(nwin):
        lo, hi = edges[i], edges[i + 1]
        # only attempt false position where the window brackets a root
        if func(lo, k) * func(hi, k) <= 0:
            r, it = false_position(func, lo, hi, k)
            if r is not None:
                # avoid recording the same root twice at shared window edges
                if not any(abs(r - rr) < 1e-6 for rr in roots):
                    roots.append(r)
                    iters.append(it)
    return roots, iters

k = 0.15

# (1) False position applied to the WHOLE interval [0, 600] at once
whole_root, whole_it = false_position(f, 0.0, 600.0, k)
print("k =", k)
print("Whole-interval false position: single root returned =", whole_root)
print("Whole-interval false position: iterations =", whole_it)
print("Whole-interval false position: f(root) =", f(whole_root, k))

# (2) Windowed scan finds all roots
wroots, witers = windowed_scan(f, 0.0, 600.0, k)
print("Windowed scan: number of roots found =", len(wroots))
for j, (r, it) in enumerate(zip(wroots, witers)):
    print("Windowed root[%d] = %.10f  (iterations = %d, f = %.3e)" % (j, r, it, f(r, k)))

# --- Check: whole interval returns only ONE of the three roots ----------
found_whole = any(abs(whole_root - r) < 1e-6 for r in wroots)
print("Check: whole-interval root count = 1, windowed root count =", len(wroots))
print("Check: whole-interval root is one of the windowed roots =", found_whole)
# Explanation (one sentence):
# Because false position only guarantees a single bracketed crossing per call, running it
# on all of [0,600] collapses the three sign changes into just one returned root, while the
# windowed scan isolates each sign change separately and thus recovers all three roots.

# --- Plot -------------------------------------------------------------
X = np.linspace(0, 600, 1000)
plt.figure(figsize=(8, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X, k), label="f(X, k=%.2f)" % k)
plt.plot(wroots, [f(r, k) for r in wroots], "ro", label="windowed roots (all 3)")
plt.plot([whole_root], [f(whole_root, k)], "b*", ms=15, label="whole-interval root (only 1)")
plt.xlabel("X"); plt.ylabel("f(X, k)")
plt.title("Self-activating gene: false-position roots at k=0.15")
plt.legend()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2E.4.1_s5.png")
