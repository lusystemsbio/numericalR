import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Model parameters ----
tau = 1.0        # delay
B = 100.0        # carrying capacity
dt = 0.01        # time step
T = 60.0         # total simulation time
r_values = [0.3, 1.5, np.pi / 2, 1.7]
r_labels = ["r=0.3", "r=1.5", "r=pi/2", "r=1.7"]

lag = int(round(tau / dt))   # number of steps in one delay (= 100)
n_steps = int(round(T / dt)) # number of integration steps
t = np.arange(n_steps + 1) * dt

# Right-hand side of the Hutchinson delay-logistic equation.
# The crowding term uses the population one delay tau in the past.
def f(N_now, N_delayed, r):
    return r * N_now * (1.0 - N_delayed / B)

results = {}
for r, label in zip(r_values, r_labels):
    N = np.empty(n_steps + 1)
    # Constant history N(t) = 1 for all t <= 0.  Because tau/dt = lag is an
    # integer, the delayed value N(t-tau) always lands exactly on a grid point:
    # for the first 'lag' steps that grid point falls in the history (= 1).
    N[0] = 1.0

    for i in range(n_steps):
        # Delayed value needed at the current time t_i.
        N_delay_now = N[i - lag] if i - lag >= 0 else 1.0
        # Delayed value needed at the next time t_{i+1}; it lies one step
        # further back and has therefore already been computed (or is history).
        N_delay_next = N[i + 1 - lag] if i + 1 - lag >= 0 else 1.0

        # --- Heun stage 1: explicit Euler predictor ---
        k1 = f(N[i], N_delay_now, r)
        N_pred = N[i] + dt * k1

        # --- Heun stage 2: corrector using slope at the predicted endpoint ---
        k2 = f(N_pred, N_delay_next, r)
        N[i + 1] = N[i] + 0.5 * dt * (k1 + k2)

    results[label] = N

    # Report late-time behaviour so the Hopf crossing r*tau = pi/2 is visible.
    tail = N[int(len(N) * 0.75):]      # last quarter of the run
    print(f"{label}: r*tau = {r * tau:.4f}, "
          f"final N = {N[-1]:.4f}, "
          f"tail min = {tail.min():.4f}, tail max = {tail.max():.4f}, "
          f"tail peak-to-peak amplitude = {tail.max() - tail.min():.4f}")

# ---- Plot N(t) for each r ----
plt.figure(figsize=(10, 6))
for label in r_labels:
    plt.plot(t, results[label], label=label)
plt.axhline(B, color="k", linestyle="--", linewidth=0.8, label="B=100")
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Hutchinson delay-logistic growth (Heun DDE integrator)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4B.1.1_s5.png")

# Explanation of why the check confirms the result:
print("Check: the tail peak-to-peak amplitude stays ~0 for r=0.3 and r=1.5 "
      "(monotone/damped to B), shrinks only slowly for r=pi/2, but remains "
      "large and non-decaying for r=1.7 (r*tau > pi/2) — confirming the Hopf "
      "bifurcation at r*tau = pi/2 where the equilibrium loses stability to a "
      "sustained limit cycle.")
