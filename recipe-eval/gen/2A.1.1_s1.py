import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Model parameters
g = 50.0    # transcription rate, nM/min
k = 0.1     # degradation rate, per min

# Steady state: set dX/dt = 0  ->  g - k*Xss = 0  ->  Xss = g/k
Xss = g / k
print("Steady state g/k =", Xss, "nM")

# Time grid (min) and initial conditions (nM)
t = np.linspace(0, 80, 801)
X0_list = [300, 400, 500, 600, 700]

# Prepare the plot
plt.figure(figsize=(8, 5))

# Evaluate the exact solution explicitly for each initial condition
for X0 in X0_list:
    # X(t) = g/k + (X0 - g/k) * exp(-k t): steady value plus a decaying offset
    offset = X0 - Xss          # signed distance from steady state at t = 0
    decay = np.exp(-k * t)     # exponential relaxation factor
    X = Xss + offset * decay   # direct evaluation, no numerical integration
    plt.plot(t, X, label="X0 = %d nM" % X0)

    # ---- check for this curve ----
    starts_at_X0 = np.isclose(X[0], X0)               # begins at its own X0
    approaches_ss = np.isclose(X[-1], Xss, atol=1.0)  # ends near steady state
    direction = "decays to" if X0 > Xss else ("rises to" if X0 < Xss else "stays at")
    print("X0 = %3d nM: X(0) = %.4f (starts at X0: %s), "
          "X(80) = %.4f (approaches g/k: %s), curve %s steady state"
          % (X0, X[0], starts_at_X0, X[-1], approaches_ss, direction))

# Mark the steady state
plt.axhline(Xss, color="k", linestyle="--", linewidth=1, label="steady state g/k = %d nM" % Xss)

plt.xlabel("time (min)")
plt.ylabel("X (nM)")
plt.title("Constitutive gene expression: dX/dt = g - k*X")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2A.1.1_s1.png")

# Explanation of why the check confirms the result:
print("Check confirms the result because each curve reproducing its own X0 at t=0 "
      "and converging to g/k=500 nM (from above if X0>500, from below if X0<500) is "
      "exactly the behavior the exact solution X(t)=g/k+(X0-g/k)exp(-k t) predicts.")
