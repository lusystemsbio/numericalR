import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model parameters ---
tau = 1.0          # delay time
dt = 0.01          # Euler time step
t_end = 40.0       # final integration time
N_history = 1.0    # constant history value N(t) = 1 for t <= 0
rates = [-0.3, -1.4, -1.7]   # growth rates to test

# Number of steps making up one delay interval (delayed value is tau/dt steps back)
delay_steps = int(round(tau / dt))
# Number of forward integration steps from t = 0 to t = t_end
n_steps = int(round(t_end / dt))

def integrate_dde(r):
    # Full trajectory INCLUDING the history interval [-tau, 0].
    # Indices 0 .. delay_steps-1 hold the history, index delay_steps corresponds to t = 0.
    total_len = delay_steps + n_steps + 1
    N = np.empty(total_len)
    t = np.empty(total_len)

    # Fill the history interval with the constant history N = 1.
    for i in range(delay_steps + 1):
        N[i] = N_history
        t[i] = -tau + i * dt

    # Explicit Euler for the DDE: N_next = N + dt * r * N_delayed
    # where N_delayed is the value tau/dt steps in the past.
    idx0 = delay_steps  # array index of t = 0
    for k in range(n_steps):
        i = idx0 + k               # current index
        N_delayed = N[i - delay_steps]   # value one delay time ago
        N[i + 1] = N[i] + dt * r * N_delayed
        t[i + 1] = t[i] + dt

    return t, N

# --- Integrate for each growth rate and report a summary value ---
results = {}
for r in rates:
    t, N = integrate_dde(r)
    results[r] = (t, N)
    # Report the final value and the peak magnitude as simple diagnostics.
    final_val = N[-1]
    max_abs = np.max(np.abs(N))
    min_val = np.min(N)
    print(f"r = {r}: N(t=40) = {final_val:.6e}")
    print(f"r = {r}: max|N| over trajectory = {max_abs:.6e}")
    print(f"r = {r}: min N over trajectory = {min_val:.6e}")

# --- Report the critical rate ---
r_crit = -np.pi / 2.0
print(f"Critical growth rate -pi/2 = {r_crit:.6f}")

# --- Plot ---
plt.figure(figsize=(9, 6))
for r in rates:
    t, N = results[r]
    plt.plot(t, N, label=f"r = {r}")
plt.axhline(0.0, color="k", lw=0.5)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed exponential growth: dN/dt = r*N(t - tau), tau = 1")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.1.1_s3.png")

# Explanation of why the check confirms the result:
print("Explanation: Because r = -0.3 decays monotonically, r = -1.4 (between -pi/2 and 0) "
      "decays with damped oscillations, and r = -1.7 (past -pi/2) oscillates with growing "
      "amplitude, the trajectories bracket the critical value -pi/2 and show the qualitative "
      "change in stability of the N = 0 state exactly where theory predicts it.")
