import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model definition: self-inhibiting gene ---
# f(X) = basal transcription + repressive Hill term - linear degradation
g0 = 10.0     # basal transcription rate
g1 = 60.0     # max repressible transcription rate
Xth = 200.0   # repression threshold (nM)
n = 4         # Hill coefficient
k = 0.1       # linear degradation rate

def f(X):
    return g0 + g1 / (1.0 + (X / Xth) ** n) - k * X

# --- Step 1: locate the steady state (f(X) = 0) near 250 nM ---
# Explicit bisection rather than a black-box root finder.
a, b = 100.0, 400.0   # bracket known to contain the root (f changes sign)
fa, fb = f(a), f(b)
assert fa * fb < 0, "root not bracketed"
for _ in range(200):          # iterate until interval is tiny
    m = 0.5 * (a + b)         # midpoint
    fm = f(m)
    if fa * fm <= 0:          # root is in [a, m]
        b, fb = m, fm
    else:                     # root is in [m, b]
        a, fa = m, fm
    if (b - a) < 1e-9:
        break
Xss = 0.5 * (a + b)           # steady-state estimate

# --- Step 2: slope df/dX via central finite difference ---
# df/dX ~ (f(X+h) - f(X-h)) / (2h)
h = 1e-4
dfdX = (f(Xss + h) - f(Xss - h)) / (2.0 * h)

# --- Step 3: stability verdict ---
# Negative slope => perturbations decay => stable.
sign = "negative" if dfdX < 0 else ("positive" if dfdX > 0 else "zero")
stable = dfdX < 0
verdict = "STABLE" if stable else "UNSTABLE"

# --- Separate confirmation check ---
# Confirm df/dX at the steady state (~250 nM) is negative.
confirmed = dfdX < 0

print(f"Steady-state value Xss = {Xss:.6f} nM")
print(f"f(Xss) (residual, should be ~0) = {f(Xss):.3e}")
print(f"df/dX at steady state = {dfdX:.6f}")
print(f"Sign of df/dX = {sign}")
print(f"Stability verdict = {verdict}")
print(f"Confirmation check: df/dX < 0 is {confirmed} -> steady state is stable")

# Why this check confirms the result: a negative slope of f at the fixed point
# means a small increase in X makes dX/dt negative (and a small decrease makes it
# positive), so the system is pushed back toward the steady state, i.e. it is stable.

# --- Plot for visualization ---
X = np.linspace(0, 500, 1000)
plt.figure(figsize=(8, 5))
plt.axhline(0, color="gray", lw=0.8)
plt.plot(X, f(X), label="f(X)")
plt.plot(Xss, 0.0, "ro", label=f"steady state {Xss:.1f} nM ({verdict})")
plt.xlabel("X (nM)")
plt.ylabel("dX/dt = f(X)")
plt.title("Self-inhibiting gene: steady state and stability")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.1.1_s5.png")
