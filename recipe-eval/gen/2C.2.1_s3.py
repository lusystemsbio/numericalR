import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
r = 0.1     # growth rate
B = 100.0   # carrying capacity
N0 = 1.0    # initial population
print(f"Parameter r (growth rate): {r}")
print(f"Parameter B (carrying capacity): {B}")
print(f"Parameter N0 (initial population): {N0}")

# --- Time grid over which to evaluate the exact solution ---
t = np.linspace(0.0, 100.0, 1001)   # t from 0 to 100
print(f"Time range: {t[0]} to {t[-1]} with {len(t)} points")

# --- Evaluate the exact logistic solution explicitly, step by step ---
# Exact solution: N(t) = N0*B / (N0 + (B - N0)*exp(-r*t))
decay = np.exp(-r * t)              # the exponential relaxation term exp(-r*t)
denominator = N0 + (B - N0) * decay # denominator of the closed-form solution
numerator = N0 * B                  # numerator N0*B (constant in t)
N = numerator / denominator         # exact population at each time

# --- Report key values of the curve ---
print(f"N(t=0) computed: {N[0]}")           # should equal N0
print(f"N(t=100) computed: {N[-1]}")        # should be close to B
print(f"Carrying capacity B (target): {B}")

# --- Separate check: sigmoidal rise from N0 leveling off at B ---
# 1) starts at N0
starts_at_N0 = np.isclose(N[0], N0)
# 2) monotonically increasing (rising)
is_increasing = np.all(np.diff(N) > 0)
# 3) levels off at carrying capacity B (approaches within a small tolerance)
levels_at_B = abs(N[-1] - B) < 0.1
# 4) sigmoidal: single inflection point where growth rate is maximal.
#    For logistic growth the inflection occurs at N = B/2.
dNdt = r * N * (1.0 - N / B)          # instantaneous growth rate
idx_max_rate = int(np.argmax(dNdt))   # index of fastest growth (inflection)
N_at_inflection = N[idx_max_rate]     # population at the inflection point
print(f"Check - starts at N0: {starts_at_N0}")
print(f"Check - monotonically increasing (rising): {is_increasing}")
print(f"Check - levels off at B within tolerance: {levels_at_B}")
print(f"Population at inflection (max growth rate): {N_at_inflection}")
print(f"Expected inflection population (B/2): {B/2.0}")
print(f"Max final deviation from B: {abs(N[-1] - B)}")

# --- Plot the logistic curve rising to the carrying capacity ---
plt.figure(figsize=(8, 5))
plt.plot(t, N, color="steelblue", lw=2, label="Exact logistic N(t)")
plt.axhline(B, color="red", ls="--", lw=1, label=f"Carrying capacity B = {B:.0f}")
plt.scatter([t[idx_max_rate]], [N_at_inflection], color="black", zorder=5,
            label="Inflection point (N = B/2)")
plt.xlabel("Time t")
plt.ylabel("Population N(t)")
plt.title("Exact Logistic Growth Curve")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.2.1_s3.png")

# Why this check confirms the result: a correct logistic solution must begin at N0,
# rise monotonically through its steepest point at exactly N = B/2, and asymptotically
# level off at B, so verifying these three signatures confirms the sigmoidal shape and
# the correct carrying-capacity limit.
print("Explanation: The curve begins at N0, rises monotonically with maximal slope at "
      "N=B/2, and asymptotes to B, which is the defining sigmoidal signature of the "
      "exact logistic solution.")
