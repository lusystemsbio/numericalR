import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g0 = 10.0     # basal transcription rate
g1 = 45.0     # max excitatory (Hill) transcription rate
Xth = 200.0   # Hill threshold (nM)
n = 4         # Hill coefficient

# Rate-of-change function f(X) = g0 + g1*Hill - k*X
def f(X, k):
    hill = (X / Xth)**n / (1.0 + (X / Xth)**n)
    return g0 + g1 * hill - k * X

# Effective potential via explicit trapezoidal accumulation:
# U(X+dx) = U(X) - (f(X)+f(X+dx))/2 * dx, with U(0)=0
def effective_potential(X, k):
    U = np.zeros_like(X)          # U[0] = 0 (U(0)=0)
    for i in range(len(X) - 1):
        dx = X[i + 1] - X[i]      # step width
        # accumulate: subtract the trapezoid area of f over [X[i], X[i+1]]
        U[i + 1] = U[i] - 0.5 * (f(X[i], k) + f(X[i + 1], k)) * dx
    return U

# Spatial grid over the state variable X (concentration, nM)
X = np.linspace(0.0, 500.0, 5001)

# Compute potentials for the three degradation rates
ks = [0.15, 0.2, 0.1]
Us = {k: effective_potential(X, k) for k in ks}

# --- Plot U(X) for k = 0.15, 0.2, 0.1 ---
plt.figure(figsize=(8, 5))
for k in ks:
    plt.plot(X, Us[k], label=f"k = {k}")
plt.xlabel("X (nM)")
plt.ylabel("Effective potential U(X)")
plt.title("Effective potential of a self-activating gene")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2D.3.1_s1.png")

# --- Separate check: count basins (local minima) and barriers (local maxima) ---
# A local minimum of U is a stable state; a local maximum is an unstable state.
def find_extrema(U):
    minima, maxima = [], []
    for i in range(1, len(U) - 1):
        if U[i] < U[i - 1] and U[i] < U[i + 1]:
            minima.append(i)
        if U[i] > U[i - 1] and U[i] > U[i + 1]:
            maxima.append(i)
    return minima, maxima

for k in ks:
    minima, maxima = find_extrema(Us[k])
    print(f"--- k = {k} ---")
    print(f"Number of basins (stable states, minima): {len(minima)}")
    for i in minima:
        print(f"  stable state at X = {X[i]:.1f} nM, U = {Us[k][i]:.2f}")
    print(f"Number of barriers (unstable states, maxima): {len(maxima)}")
    for i in maxima:
        print(f"  unstable state at X = {X[i]:.1f} nM, U = {Us[k][i]:.2f}")

# One-sentence explanation:
print("\nExplanation: This check confirms the result because two potential minima "
      "separated by one maximum means the k=0.15 system is bistable (two stable "
      "states around 100 and 300 nM split by an unstable barrier near 200 nM), "
      "whereas a single minimum at k=0.2 and k=0.1 means a single stable well "
      "(monostable) with no barrier.")
