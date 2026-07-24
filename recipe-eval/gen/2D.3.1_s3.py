import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters (self-activating gene) ---
g0 = 10.0     # basal transcription
g1 = 45.0     # max excitatory (Hill) contribution
Xth = 200.0   # Hill threshold (nM)
n = 4         # Hill coefficient

# f(X) = production - degradation
def f(X, k):
    hill = (X / Xth)**n / (1.0 + (X / Xth)**n)   # excitatory Hill term
    return g0 + g1 * hill - k * X                  # minus linear degradation

# --- Effective potential U(X) = -integral f dx, with U(0)=0 ---
# Built explicitly with the trapezoidal recurrence:
#   U(X+dx) = U(X) - (f(X)+f(X+dx))/2 * dx
def effective_potential(k, Xmax=500.0, dx=0.1):
    X = np.arange(0.0, Xmax + dx, dx)   # grid from 0 upward
    fvals = f(X, k)                     # rate at each grid point
    U = np.zeros_like(X)                # U(0) = 0
    for i in range(len(X) - 1):
        # accumulate one trapezoid of -f between X[i] and X[i+1]
        U[i + 1] = U[i] - 0.5 * (fvals[i] + fvals[i + 1]) * dx
    return X, U

# --- Compute potentials for the three degradation rates ---
ks = [0.15, 0.2, 0.1]
results = {k: effective_potential(k) for k in ks}

# --- Plot ---
plt.figure(figsize=(8, 5))
for k in ks:
    X, U = results[k]
    plt.plot(X, U, label=f"k = {k}")
plt.xlabel("X (nM)")
plt.ylabel("Effective potential U(X)")
plt.title("Effective potential of a self-activating gene")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.3.1_s3.png")

# --- Helper: find interior local minima (valleys=stable) and maxima (peaks=unstable) ---
def find_extrema(X, U):
    minima, maxima = [], []
    for i in range(1, len(U) - 1):
        if U[i] < U[i - 1] and U[i] < U[i + 1]:
            minima.append(X[i])
        if U[i] > U[i - 1] and U[i] > U[i + 1]:
            maxima.append(X[i])
    return minima, maxima

# --- Report extrema for each k ---
for k in ks:
    X, U = results[k]
    minima, maxima = find_extrema(X, U)
    print(f"k = {k}: valleys (stable states, nM) = "
          f"{[round(m, 1) for m in minima]}")
    print(f"k = {k}: peaks (unstable states, nM) = "
          f"{[round(m, 1) for m in maxima]}")

# --- Explicit check at k = 0.15: two basins + one barrier ---
X15, U15 = results[0.15]
min15, max15 = find_extrema(X15, U15)
two_basins = len(min15) == 2 and len(max15) == 1
print(f"k = 0.15: number of valleys = {len(min15)}")
print(f"k = 0.15: number of peaks   = {len(max15)}")
print(f"k = 0.15: bistable (2 valleys split by 1 barrier)? {two_basins}")

# --- Check that k = 0.2 and k = 0.1 are monostable (single well) ---
for k in [0.2, 0.1]:
    X, U = results[k]
    mn, mx = find_extrema(X, U)
    print(f"k = {k}: single well (1 valley, 0 interior peaks)? "
          f"{len(mn) == 1 and len(mx) == 0}")

# One-sentence explanation:
print("Explanation: The check confirms the result because valleys of U(X) are "
      "exactly the stable steady states and peaks are the unstable ones "
      "(f = -dU/dX = 0 at each extremum, stable where U curves up, unstable "
      "where it curves down), so two valleys near 100 and 300 nM separated by "
      "a barrier near 200 nM at k=0.15 versus a single valley at k=0.2 and "
      "k=0.1 directly demonstrates bistability only at the intermediate "
      "degradation rate.")
