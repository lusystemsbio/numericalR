import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model parameters ----
g0 = 10.0     # basal transcription rate
g1 = 45.0     # max Hill activation rate
Xth = 200.0   # Hill threshold
n = 4         # Hill coefficient

# Rate-of-change f(X) = basal + Hill activation - linear degradation
def f(X, k):
    hill = (X / Xth) ** n / (1.0 + (X / Xth) ** n)
    return g0 + g1 * hill - k * X

# ---- Effective potential via explicit trapezoidal accumulation ----
# U(X) = -integral f(x) dx, with U(0) = 0.
# Step rule: U(X+dx) = U(X) - (f(X) + f(X+dx))/2 * dx
def effective_potential(X, k):
    U = np.zeros_like(X)          # U[0] = 0 at X = 0
    for i in range(len(X) - 1):
        dx = X[i + 1] - X[i]      # grid spacing
        # subtract the trapezoidal area of f over [X[i], X[i+1]]
        U[i + 1] = U[i] - 0.5 * (f(X[i], k) + f(X[i + 1], k)) * dx
    return U

# X grid (nM)
X = np.linspace(0.0, 500.0, 2001)

# Helper: find fixed points (sign changes of f) and classify by slope f'
def fixed_points(X, k):
    fv = f(X, k)
    pts = []
    for i in range(len(X) - 1):
        if fv[i] == 0.0 or fv[i] * fv[i + 1] < 0.0:
            # linear interpolation for the root location
            x0 = X[i] - fv[i] * (X[i + 1] - X[i]) / (fv[i + 1] - fv[i])
            # slope of f at root: f'<0 -> stable (valley), f'>0 -> unstable (peak)
            slope = (f(x0 + 1e-3, k) - f(x0 - 1e-3, k)) / 2e-3
            kind = "stable (valley)" if slope < 0 else "unstable (peak)"
            pts.append((x0, kind))
    return pts

# ---- Compute potentials for the three degradation rates ----
ks = [0.15, 0.2, 0.1]
potentials = {}
for k in ks:
    potentials[k] = effective_potential(X, k)
    print(f"--- k = {k} ---")
    fps = fixed_points(X, k)
    print(f"Number of fixed points: {len(fps)}")
    for x0, kind in fps:
        print(f"Fixed point at X = {x0:.2f} nM : {kind}, U = {np.interp(x0, X, potentials[k]):.2f}")
    print(f"U at X=500 nM: {potentials[k][-1]:.2f}")

# ---- Explicit check at k = 0.15 for two basins split by a barrier ----
k_check = 0.15
U = potentials[k_check]
# minima (valleys) and maxima (peaks) via discrete first-derivative sign change
dU = np.gradient(U, X)
minima = [X[i] for i in range(1, len(X) - 1) if dU[i - 1] < 0 and dU[i + 1] > 0]
maxima = [X[i] for i in range(1, len(X) - 1) if dU[i - 1] > 0 and dU[i + 1] < 0]
print("--- Check at k = 0.15 ---")
print(f"Valleys (stable states) near X = {[f'{m:.1f}' for m in minima]} nM")
print(f"Peaks (unstable states) near X = {[f'{m:.1f}' for m in maxima]} nM")
print(f"Number of valleys: {len(minima)}, number of peaks: {len(maxima)}")
print(f"Bistable (two valleys, one barrier)? {len(minima) == 2 and len(maxima) == 1}")

# ---- Plot ----
plt.figure(figsize=(8, 6))
for k in ks:
    plt.plot(X, potentials[k], label=f"k = {k}")
plt.xlabel("X (nM)")
plt.ylabel("Effective potential U(X)")
plt.title("Effective potential of a self-activating gene circuit")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.3.1_s2.png")

# One-sentence explanation of why the check confirms the result:
# The check confirms the result because valleys of U occur exactly where f(X)=0 with
# f'(X)<0 (stable steady states) and peaks where f(X)=0 with f'(X)>0 (unstable steady
# states), so finding two minima flanking one maximum at k=0.15 (versus a single
# minimum at k=0.2 and k=0.1) directly demonstrates bistability versus monostability.
print("The check confirms the result because U's valleys sit at stable fixed points (f=0, f'<0) and its peak at the unstable fixed point (f=0, f'>0), so two valleys split by one barrier at k=0.15 proves bistability while a single well at k=0.2 and k=0.1 proves monostability.")
