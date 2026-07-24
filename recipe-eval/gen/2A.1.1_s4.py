import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g = 50.0    # transcription rate (nM/min)
k = 0.1     # degradation rate (per min)

# Steady state: at dX/dt = 0 -> g - k*X = 0 -> X = g/k
Xss = g / k
print(f"Steady state g/k = {Xss:.4f} nM")

# --- Time grid and initial conditions ---
t = np.linspace(0, 80, 801)   # t from 0 to 80 min
X0_list = [300.0, 400.0, 500.0, 600.0, 700.0]

# --- Direct evaluation of the exact solution (explicit, not a routine) ---
# X(t) = g/k + (X0 - g/k) * exp(-k*t)
solutions = {}
for X0 in X0_list:
    offset = X0 - Xss                 # initial displacement from steady state
    decay = np.exp(-k * t)            # exponential relaxation factor
    X = Xss + offset * decay          # combine steady state + decaying transient
    solutions[X0] = X

# --- Plot one curve per initial condition ---
plt.figure(figsize=(8, 5))
for X0 in X0_list:
    plt.plot(t, solutions[X0], label=f"X0 = {X0:.0f} nM")
plt.axhline(Xss, color="black", linestyle="--", label=f"steady state g/k = {Xss:.0f} nM")
plt.xlabel("time (min)")
plt.ylabel("X (nM)")
plt.title("Constitutive gene expression: relaxation to steady state")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2A.1.1_s4.png")

# --- Check: start value, final value, and direction of approach ---
print("\nCheck for each initial condition:")
all_ok = True
for X0 in X0_list:
    X = solutions[X0]
    start_val = X[0]                  # should equal X0
    end_val = X[-1]                   # should be near steady state
    print(f"X0 = {X0:.0f} nM: start = {start_val:.4f} nM, "
          f"end (t=80) = {end_val:.4f} nM, distance to g/k = {abs(end_val - Xss):.4f} nM")

    # verify it begins at its own X0
    starts_at_X0 = np.isclose(start_val, X0)
    # verify it approaches the steady state
    approaches_ss = abs(end_val - Xss) < abs(X0 - Xss)
    # verify correct direction: above 500 decays, below 500 rises
    if X0 > Xss:
        direction_ok = end_val < X0 and end_val > Xss
        direction = "decays toward g/k"
    elif X0 < Xss:
        direction_ok = end_val > X0 and end_val < Xss
        direction = "rises toward g/k"
    else:
        direction_ok = np.isclose(end_val, Xss)
        direction = "stays at g/k"
    print(f"    begins at X0: {starts_at_X0}, approaches g/k: {approaches_ss}, "
          f"{direction}: {direction_ok}")
    all_ok = all_ok and starts_at_X0 and approaches_ss and direction_ok

print(f"\nAll checks passed: {all_ok}")

# Explanation: the check confirms the result because the exact solution must reduce
# to X0 at t=0 (since exp(0)=1) and to g/k as t grows (since exp(-k*t)->0), so
# matching those endpoints with the right monotonic direction verifies both the
# initial condition and the guaranteed convergence to the unique steady state.
print("\nWhy this check confirms the result: at t=0 exp(-k*t)=1 forces X=X0 and as "
      "t grows exp(-k*t)->0 forces X->g/k, so verifying each curve starts at its own "
      "X0 and monotonically approaches 500 nM confirms both the initial condition and "
      "the correct relaxation to the unique steady state.")
