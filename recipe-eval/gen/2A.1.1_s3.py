import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g = 50.0    # transcription rate, nM/min
k = 0.1     # degradation rate, per min

# --- Steady state: set dX/dt = 0 => g - k*X = 0 => X = g/k ---
steady_state = g / k
print(f"Steady state g/k: {steady_state:.4f} nM")

# --- Time grid and initial conditions ---
t = np.linspace(0, 80, 400)          # time from 0 to 80 min
X0_list = [300, 400, 500, 600, 700]  # initial values in nM

# --- Direct evaluation of the exact solution (no integration) ---
# X(t) = g/k + (X0 - g/k)*exp(-k*t)
solutions = {}
for X0 in X0_list:
    X = steady_state + (X0 - steady_state) * np.exp(-k * t)  # closed-form value at each t
    solutions[X0] = X

# --- Check: each curve starts at its own X0 and approaches g/k ---
print("\nCheck of start and end values:")
for X0 in X0_list:
    X = solutions[X0]
    start = X[0]            # value at t = 0
    end = X[-1]             # value at t = 80 min (near steady state)
    direction = "decays to" if X0 > steady_state else ("rises to" if X0 < steady_state else "stays at")
    print(f"X0 = {X0:>3d} nM: starts at {start:.4f} nM, ends at {end:.4f} nM, {direction} steady state")
    # Verify programmatically
    assert np.isclose(start, X0), "start does not equal X0"
    assert abs(end - steady_state) < 1.0, "end does not approach steady state"

# --- Plot ---
plt.figure(figsize=(8, 5))
for X0 in X0_list:
    plt.plot(t, solutions[X0], label=f"X0 = {X0} nM")
plt.axhline(steady_state, color="black", linestyle="--", label=f"steady state g/k = {steady_state:.0f} nM")
plt.xlabel("time (min)")
plt.ylabel("X (nM)")
plt.title("Constitutive gene expression: relaxation to steady state")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2A.1.1_s3.png")

# The check confirms the result because every curve's value at t=0 equals its own X0
# while all curves converge to the same g/k = 500 nM, with those above 500 decaying and
# those below rising, exactly matching the sign of (X0 - g/k) in the exponential term.
print("\nThe check confirms the result because each curve begins exactly at its own X0 and "
      "all approach g/k = 500 nM, with curves above decaying and below rising, matching the "
      "sign of (X0 - g/k) in the decaying exponential.")
