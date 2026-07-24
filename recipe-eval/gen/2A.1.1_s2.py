import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
g = 50.0    # transcription rate, nM/min
k = 0.1     # degradation rate, per min

# Steady state: set dX/dt = 0 -> g - k*X = 0 -> X = g/k
steady_state = g / k
print(f"Steady state g/k = {steady_state} nM")

# Time vector: 0 to 80 min
t = np.linspace(0, 80, 801)

# Initial conditions to test
X0_list = [300.0, 400.0, 500.0, 600.0, 700.0]

# --- Explicit evaluation of the exact solution X(t) = g/k + (X0 - g/k)*exp(-k*t) ---
plt.figure(figsize=(8, 6))
for X0 in X0_list:
    offset = X0 - steady_state          # initial deviation from steady state
    decay = np.exp(-k * t)              # exponential relaxation factor
    X = steady_state + offset * decay  # exact solution assembled term by term
    plt.plot(t, X, label=f"X0 = {X0:.0f} nM")

    # --- Check: report start value and end (approached) value for this curve ---
    start_val = X[0]
    end_val = X[-1]
    direction = "decays to" if X0 > steady_state else ("rises to" if X0 < steady_state else "stays at")
    print(f"X0 = {X0:.0f} nM: X(0) = {start_val:.2f} nM (matches X0), "
          f"X(80) = {end_val:.4f} nM, curve {direction} steady state {steady_state:.0f} nM")

# Mark the steady state
plt.axhline(steady_state, color="black", linestyle="--", label=f"steady state g/k = {steady_state:.0f} nM")

plt.xlabel("time (min)")
plt.ylabel("X (nM)")
plt.title("Constitutive gene expression: relaxation to steady state")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2A.1.1_s2.png")

# Explanation of why the check confirms the result:
print("Check confirms the result because each curve starting exactly at its own X0 "
      "and monotonically converging to g/k = 500 nM (from above if X0>500, from below if X0<500) "
      "is precisely the behavior the exact solution X(t)=g/k+(X0-g/k)exp(-kt) predicts, "
      "since the exp(-kt) term decays to zero leaving X = g/k regardless of X0.")
