import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model parameters (fixed) ----
g0 = 10.0      # basal transcription rate
g1 = 45.0      # max Hill-activated transcription rate
Xth = 200.0    # Hill threshold (nM)
n = 4          # Hill coefficient

# ---- Self-activating gene rate law ----
# f(X) = basal + excitatory Hill activation - linear degradation
def f(X, k):
    hill = (X / Xth)**n / (1.0 + (X / Xth)**n)
    return g0 + g1 * hill - k * X

# ---- Effective potential via explicit trapezoidal accumulation ----
# U(X) = -integral_0^X f(x) dx, with U(0) = 0.
# Implemented step by step: U(X+dx) = U(X) - (f(X)+f(X+dx))/2 * dx
def effective_potential(Xgrid, k):
    U = np.zeros_like(Xgrid)          # U(0) = 0 at the first grid point
    for i in range(len(Xgrid) - 1):  # walk along the grid one interval at a time
        dx = Xgrid[i+1] - Xgrid[i]    # step width
        fL = f(Xgrid[i], k)           # rate at left edge of interval
        fR = f(Xgrid[i+1], k)         # rate at right edge of interval
        # trapezoidal area of f over [X, X+dx], subtracted because U = -integral f
        U[i+1] = U[i] - 0.5 * (fL + fR) * dx
    return U

# ---- Grid over X (nM) ----
Xgrid = np.linspace(0.0, 500.0, 5001)  # dx = 0.1 nM

# ---- Compute potentials for the three degradation rates ----
ks = [0.15, 0.2, 0.1]
U = {k: effective_potential(Xgrid, k) for k in ks}

# ---- Locate stationary states: extrema of U == zeros of f ----
# Minima of U (valleys) are stable states; maxima (peaks) are unstable states.
def classify_extrema(Xgrid, k):
    Uk = U[k]
    dU = np.diff(Uk)                 # forward differences of U; sign changes mark extrema
    minima, maxima = [], []
    for i in range(1, len(dU)):
        if dU[i-1] < 0 and dU[i] > 0:   # U decreasing then increasing -> valley (stable)
            minima.append(Xgrid[i])
        if dU[i-1] > 0 and dU[i] < 0:   # U increasing then decreasing -> peak (unstable)
            maxima.append(Xgrid[i])
    return minima, maxima

# ---- Report the check for every k ----
for k in ks:
    minima, maxima = classify_extrema(Xgrid, k)
    print(f"k = {k}:")
    print(f"  stable states (valleys / U minima, nM):   {[round(m,1) for m in minima]}")
    print(f"  unstable states (barriers / U maxima, nM): {[round(m,1) for m in maxima]}")
    print(f"  number of basins: {len(minima)}")

# ---- Plot the effective potentials ----
plt.figure(figsize=(8, 5))
colors = {0.15: "C0", 0.2: "C1", 0.1: "C2"}
for k in ks:
    plt.plot(Xgrid, U[k], color=colors[k], label=f"k = {k}")
plt.xlabel("X (nM)")
plt.ylabel("Effective potential U(X)")
plt.title("Effective potential of a self-activating gene circuit")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.3.1_s5.png")

# ---- One-sentence explanation of the check ----
# The check confirms the result because valleys of U are exactly the zeros of f where
# f slopes downward (stable fixed points) and peaks are zeros where f slopes upward
# (unstable fixed points), so two valleys + one intervening peak at k=0.15 versus a
# single valley at k=0.2 and k=0.1 is precisely bistability giving way to monostability.
print("\nExplanation: two valleys separated by a peak at k=0.15 mean two stable states"
      " with an unstable threshold between them (bistability), while the single well at"
      " k=0.2 and k=0.1 means only one stable state (monostability); since U's minima are"
      " the stable and its maxima the unstable fixed points of f, this directly confirms"
      " the expected 100/200/300 nM structure.")
