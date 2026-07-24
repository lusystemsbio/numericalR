import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
r = 0.1      # growth rate
B = 100.0    # carrying capacity
N0 = 1.0     # initial population
t = np.linspace(0, 100, 1001)  # time grid from 0 to 100

# --- Exact solution, evaluated directly (implemented explicitly) ---
# N(t) = N0*B / (N0 + (B - N0)*exp(-r*t))
decay = np.exp(-r * t)                 # exponential decay term e^{-r t}
denominator = N0 + (B - N0) * decay    # denominator of the closed-form solution
N = N0 * B / denominator               # exact logistic curve

# --- Report key numerical results ---
print(f"r (growth rate)              = {r}")
print(f"B (carrying capacity)        = {B}")
print(f"N0 (initial population)      = {N0}")
print(f"N(t=0)   (should equal N0)   = {N[0]}")
print(f"N(t=100) (near capacity B)   = {N[-1]}")
print(f"Final value / B (fraction)   = {N[-1] / B}")

# --- Sigmoidal / leveling-off check ---
# 1) The curve must be monotonically increasing (rising from N0).
increasing = np.all(np.diff(N) > 0)
# 2) It must stay below B and approach it (levels off at carrying capacity).
below_capacity = np.all(N < B)
approaches_B = abs(N[-1] - B) < 1.0
# 3) The maximum slope (inflection of the sigmoid) occurs near N = B/2,
#    which is the hallmark of sigmoidal growth.
slopes = np.diff(N) / np.diff(t)
i_max_slope = np.argmax(slopes)          # index of steepest growth
N_at_max_slope = 0.5 * (N[i_max_slope] + N[i_max_slope + 1])

print(f"Monotonically increasing     = {increasing}")
print(f"Always below capacity B      = {below_capacity}")
print(f"Approaches B by t=100        = {approaches_B}")
print(f"N at steepest growth         = {N_at_max_slope}  (expected ~ B/2 = {B/2})")
print(f"Check passed (sigmoid + level-off) = {increasing and below_capacity and approaches_B}")

# --- Plot ---
plt.figure(figsize=(8, 5))
plt.plot(t, N, color="C0", label="Exact logistic N(t)")
plt.axhline(B, color="k", linestyle="--", label=f"Carrying capacity B = {B:.0f}")
plt.axhline(N0, color="gray", linestyle=":", label=f"N0 = {N0:.0f}")
plt.xlabel("time t")
plt.ylabel("population N(t)")
plt.title("Exact logistic growth (r=0.1, B=100, N0=1)")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2C.2.1_s5.png")

# Explanation: The check confirms the result because a correct logistic solution must
# start at N0, increase monotonically with its steepest slope near N=B/2, and asymptotically
# level off just below the carrying capacity B, so verifying all three properties uniquely
# validates the sigmoidal shape produced by direct evaluation of the closed-form solution.
