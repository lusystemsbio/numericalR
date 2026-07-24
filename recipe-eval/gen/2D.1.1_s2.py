import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model definition ----
# Self-inhibiting gene: basal transcription + repressive Hill term - linear degradation
def f(X, g0=10.0, g1=60.0, Xth=200.0, n=4.0, k=0.1):
    return g0 + g1 / (1.0 + (X / Xth) ** n) - k * X

# ---- Step 1: locate the steady state where f(X) = 0, near 250 nM ----
# Use bisection explicitly (rather than a one-step solver) to bracket the root.
lo, hi = 100.0, 400.0          # bracket that surrounds the expected root ~250
assert f(lo) > 0 and f(hi) < 0  # f is positive at low X, negative at high X
for _ in range(200):            # iterate until the bracket is tiny
    mid = 0.5 * (lo + hi)
    if f(mid) > 0:
        lo = mid                # root is to the right
    else:
        hi = mid                # root is to the left
Xss = 0.5 * (lo + hi)           # steady-state estimate

# ---- Step 2: central finite-difference estimate of df/dX at the steady state ----
h = 1e-4                                    # small step for the difference
dfdX = (f(Xss + h) - f(Xss - h)) / (2.0 * h)  # central difference formula

# ---- Step 3: linear-stability verdict from the sign of the slope ----
slope_sign = "negative" if dfdX < 0 else ("positive" if dfdX > 0 else "zero")
stable = dfdX < 0
verdict = "STABLE" if stable else "UNSTABLE"

# ---- Separate confirmation check ----
# A negative slope means a small perturbation decays back to the fixed point,
# which is exactly the condition for linear stability, confirming the verdict.
confirm_negative = dfdX < 0

# ---- Print all numerical results ----
print(f"Steady-state value Xss (nM): {Xss:.6f}")
print(f"Residual f(Xss): {f(Xss):.3e}")
print(f"df/dX at steady state (central difference): {dfdX:.6e}")
print(f"Sign of df/dX: {slope_sign}")
print(f"Stability verdict: {verdict}")
print(f"Confirmation check (df/dX < 0 => stable): {confirm_negative}")
print("Why: df/dX < 0 means small deviations from the fixed point decay, so the state is stable.")

# ---- Plot f(X) with the steady state marked ----
X = np.linspace(0, 500, 600)
plt.figure(figsize=(7, 4.5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X), label="f(X) = g0 + g1/(1+(X/Xth)^n) - k*X")
plt.plot(Xss, 0.0, "ro", label=f"steady state ~{Xss:.1f} nM ({verdict})")
plt.xlabel("X (nM)")
plt.ylabel("f(X)")
plt.title(f"Self-inhibiting gene: df/dX = {dfdX:.3e} ({slope_sign}) -> {verdict}")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.1.1_s2.png")
