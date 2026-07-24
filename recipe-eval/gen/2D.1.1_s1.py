import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g0 = 10.0     # basal transcription rate
g1 = 60.0     # max repressible transcription rate
Xth = 200.0   # Hill threshold (nM)
n = 4.0       # Hill coefficient
k = 0.1       # linear degradation rate

# --- Rate-of-change function for the self-inhibiting gene ---
# f(X) = basal + repressive Hill term - linear degradation
def f(X):
    return g0 + g1 / (1.0 + (X / Xth)**n) - k * X

# --- Step 1: find the steady state f(X) = 0 near 250 nM ---
# Use a simple bisection on an interval bracketing the root.
a, b = 100.0, 400.0     # bracket known to contain the single root
fa = f(a)
for _ in range(200):    # iterate until the interval is tiny
    m = 0.5 * (a + b)   # midpoint
    fm = f(m)
    if fa * fm <= 0.0:  # root is in [a, m]
        b = m
    else:               # root is in [m, b]
        a = m
        fa = fm
Xss = 0.5 * (a + b)     # steady-state estimate

# --- Step 2: estimate df/dX at the steady state via central finite difference ---
h = 1e-3                                  # small step for the difference
dfdX = (f(Xss + h) - f(Xss - h)) / (2.0 * h)  # central difference slope

# --- Step 3: linear stability verdict ---
# Negative slope means perturbations decay -> stable; positive means they grow -> unstable.
slope_sign = "negative" if dfdX < 0 else ("positive" if dfdX > 0 else "zero")
stable = dfdX < 0
verdict = "STABLE" if stable else "UNSTABLE"

# --- Separate confirmation check ---
# Explicitly confirm the required condition df/dX < 0 at the ~250 nM steady state.
confirm_negative = dfdX < 0

# --- Print all numerical results ---
print(f"Steady-state value Xss (nM): {Xss:.6f}")
print(f"Residual f(Xss) (should be ~0): {f(Xss):.6e}")
print(f"Finite-difference step h: {h}")
print(f"df/dX at steady state: {dfdX:.8f}")
print(f"Sign of df/dX: {slope_sign}")
print(f"Stability verdict: {verdict}")
print(f"Confirmation check df/dX < 0 is True: {confirm_negative}")

# Explanation (one sentence):
# The check confirms the result because a negative df/dX at f(X)=0 means any small
# perturbation is pushed back toward the steady state, which is the definition of linear stability.
print("Explanation: df/dX < 0 at the root means small perturbations decay back to Xss, "
      "so the steady state is linearly stable.")

# --- Plot f(X) and mark the steady state ---
Xgrid = np.linspace(0, 500, 1000)
plt.figure(figsize=(8, 5))
plt.axhline(0, color="gray", linewidth=0.8)
plt.plot(Xgrid, f(Xgrid), label="f(X)")
plt.plot(Xss, 0.0, "ro", label=f"steady state ~ {Xss:.1f} nM ({verdict})")
plt.xlabel("X (nM)")
plt.ylabel("dX/dt = f(X)")
plt.title("Self-inhibiting gene: rate function and steady state")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.1.1_s1.png")
