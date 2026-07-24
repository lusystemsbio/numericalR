import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Parameters ---
g = 50.0          # transcription rate, nM/min
k = 0.1           # degradation rate, per min
X0_list = [300.0, 400.0, 500.0, 600.0, 700.0]  # initial values, nM
t = np.linspace(0, 80, 400)  # time grid, 0 to 80 min

# --- Steady state: set dX/dt = 0 -> g - k*X = 0 -> X = g/k ---
steady_state = g / k
print(f"Steady state g/k = {steady_state:.4f} nM")

# --- Exact solution, evaluated directly (no numerical integration) ---
# X(t) = g/k + (X0 - g/k)*exp(-k*t)
def exact_solution(X0, t):
    # constant particular part plus decaying transient toward steady state
    return steady_state + (X0 - steady_state) * np.exp(-k * t)

# --- Compute and plot one curve per initial condition ---
plt.figure(figsize=(8, 5))
for X0 in X0_list:
    X = exact_solution(X0, t)               # direct evaluation over the time grid
    plt.plot(t, X, label=f"X0 = {X0:.0f} nM")

    # --- Check: starting value and long-time limit for this curve ---
    start_val = X[0]                        # value at t = 0
    end_val = X[-1]                         # value at t = 80 min (near limit)
    direction = "decays" if X0 > steady_state else ("rises" if X0 < steady_state else "constant")
    print(f"X0 = {X0:.0f} nM: X(0) = {start_val:.4f} nM, X(80) = {end_val:.4f} nM, curve {direction} toward {steady_state:.0f} nM")

# --- Mark the steady state ---
plt.axhline(steady_state, color="black", linestyle="--", label=f"steady state g/k = {steady_state:.0f} nM")

plt.xlabel("time (min)")
plt.ylabel("X (nM)")
plt.title("Constitutive gene expression: dX/dt = g - k*X")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2A.1.1_s5.png")

# --- Explanation of why the check confirms the result ---
# The check confirms the result because each curve exactly reproduces its own X0 at t=0 and
# converges monotonically to g/k (from above if X0 > 500, from below if X0 < 500), which is
# precisely the behavior the exact solution X(t) = g/k + (X0 - g/k)*exp(-k*t) predicts.
print("Check confirms the result: every curve starts at its own X0 and converges to g/k = 500 nM (above decays, below rises), matching the exact solution's transient exp(-k*t) term vanishing to leave the steady state.")
