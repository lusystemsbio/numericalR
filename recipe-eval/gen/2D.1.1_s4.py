import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters for the self-inhibiting gene
g0 = 10.0    # basal transcription rate
g1 = 60.0    # max repressible transcription rate
Xth = 200.0  # Hill threshold (nM)
n = 4        # Hill coefficient
k = 0.1      # linear degradation rate

# Rate-of-change function: production (basal + repressive Hill) minus degradation
def f(X):
    return g0 + g1 / (1.0 + (X / Xth) ** n) - k * X

# --- Find the steady state (f(X) = 0) near 250 nM via bisection ---
# Bracket the root; f is decreasing overall so pick a sign-changing interval.
a, b = 1.0, 1000.0
assert f(a) > 0 and f(b) < 0, "Root not bracketed"
for _ in range(200):          # bisection: repeatedly halve the interval
    m = 0.5 * (a + b)
    if f(a) * f(m) <= 0:
        b = m
    else:
        a = m
Xss = 0.5 * (a + b)           # steady-state estimate

# --- Central finite-difference estimate of df/dX at the steady state ---
h = 1e-3 * Xss                # small step relative to Xss
dfdX = (f(Xss + h) - f(Xss - h)) / (2.0 * h)  # central difference

# --- Stability verdict from linear stability theory ---
stable = dfdX < 0
verdict = "STABLE" if stable else "UNSTABLE"

# Print results, each labeled on its own line
print(f"Steady-state value Xss (nM): {Xss:.6f}")
print(f"f(Xss) (should be ~0): {f(Xss):.6e}")
print(f"df/dX at Xss (central difference): {dfdX:.6e}")
print(f"Sign of df/dX: {'negative' if dfdX < 0 else 'positive'}")
print(f"Stability verdict: {verdict}")

# Separate check: confirm df/dX < 0 at the steady state (~250 nM)
check_negative = dfdX < 0
print(f"Check df/dX < 0 at ~250 nM: {check_negative}")
# This check confirms the result because a negative slope means any small
# perturbation from the steady state produces a restoring change in X that
# pushes it back, which is exactly the definition of a stable fixed point.

# --- Plot f(X) and mark the steady state ---
X = np.linspace(0, 800, 800)
plt.figure(figsize=(8, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X), label="f(X) = g0 + g1/(1+(X/Xth)^n) - k*X")
plt.plot(Xss, 0, "ro", label=f"steady state ≈ {Xss:.1f} nM ({verdict})")
plt.xlabel("X (nM)")
plt.ylabel("f(X)  (rate of change)")
plt.title("Self-inhibiting gene: steady state and stability")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.1.1_s4.png")
