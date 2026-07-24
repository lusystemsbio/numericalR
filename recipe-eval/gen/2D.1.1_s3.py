import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

# --- Model parameters ---
g0 = 10.0     # basal transcription rate
g1 = 60.0     # max repressible transcription rate
Xth = 200.0   # Hill threshold
n = 4         # Hill coefficient
k = 0.1       # linear degradation rate

# --- Rate-of-change function f(X): basal + repressive Hill - degradation ---
def f(X):
    return g0 + g1 / (1.0 + (X / Xth) ** n) - k * X

# --- Find the steady state (f(X) = 0) near 250 nM ---
# Bracket a root around the expected value and solve for f(X) = 0.
X_ss = brentq(f, 100.0, 400.0)
print(f"Steady-state value X_ss = {X_ss:.6f} nM")
print(f"Residual f(X_ss) (should be ~0) = {f(X_ss):.3e}")

# --- Central finite-difference estimate of df/dX at the steady state ---
# df/dX ~= (f(X+h) - f(X-h)) / (2h), a symmetric approximation of the slope.
h = 1e-4 * X_ss                      # small step scaled to X
dfdX = (f(X_ss + h) - f(X_ss - h)) / (2.0 * h)
print(f"Central-difference df/dX at X_ss = {dfdX:.6f}")

# --- Sign of the slope and stability verdict ---
sign = "negative" if dfdX < 0 else ("positive" if dfdX > 0 else "zero")
print(f"Sign of df/dX = {sign}")
verdict = "STABLE" if dfdX < 0 else "UNSTABLE"
print(f"Stability verdict: steady state is {verdict}")

# --- Separate check near 250 nM: confirm df/dX is negative there ---
X_check = 250.0
h2 = 1e-4 * X_check
dfdX_check = (f(X_check + h2) - f(X_check - h2)) / (2.0 * h2)
print(f"Check point X = {X_check:.1f} nM, df/dX = {dfdX_check:.6f}")
print(f"df/dX at ~250 nM is {'negative (confirms STABLE)' if dfdX_check < 0 else 'not negative'}")

# Explanation (one sentence): A negative slope means any small perturbation from
# the steady state produces a rate of change that pushes X back toward it, so the
# fixed point is stable — which is exactly what the check confirms.

# --- Plot f(X) with the steady state marked ---
X = np.linspace(0, 500, 500)
plt.figure(figsize=(7, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X), label="f(X) = g0 + g1/(1+(X/Xth)^n) - kX")
plt.plot(X_ss, 0, "ro", label=f"steady state ~ {X_ss:.1f} nM ({verdict})")
plt.xlabel("X (nM)")
plt.ylabel("dX/dt = f(X)")
plt.title("Self-inhibiting gene: steady state and stability")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.1.1_s3.png")
