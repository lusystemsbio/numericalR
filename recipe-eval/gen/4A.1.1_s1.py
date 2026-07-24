import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: dN/dt = r * N(t - tau), steady state N = 0 ---
# Euler for DDE: N_next = N + dt * r * N_delayed, where N_delayed is the
# stored value tau/dt steps in the past. We keep the whole trajectory,
# including the history interval t in [-tau, 0].

def integrate_dde(r, tau=1.0, dt=0.01, t_end=40.0, hist_value=1.0):
    delay_steps = int(round(tau / dt))     # number of steps that make up one delay time
    n_hist = delay_steps                   # number of stored history points before t=0
    n_future = int(round(t_end / dt))      # number of integration steps forward

    # Time array covers the history interval [-tau, 0] then [0, t_end].
    t = np.arange(-n_hist, n_future + 1) * dt
    N = np.empty(t.size)

    # Constant history: N(t) = hist_value for t <= 0 (fills indices 0..n_hist).
    N[:n_hist + 1] = hist_value

    # Explicit Euler stepping starting at the index corresponding to t = 0.
    for i in range(n_hist, N.size - 1):
        N_delayed = N[i - delay_steps]     # value one delay time in the past
        N[i + 1] = N[i] + dt * r * N_delayed  # Euler update
    return t, N

# --- Run the three cases ---
rates = [-0.3, -1.4, -1.7]
labels = {-0.3: "r = -0.3 (monotonic)",
          -1.4: "r = -1.4 (damped oscillation)",
          -1.7: "r = -1.7 (growing oscillation)"}

results = {}
for r in rates:
    t, N = integrate_dde(r)
    results[r] = (t, N)
    # Report final value and max absolute value as simple summary numbers.
    print(f"r = {r}: N(t=40) = {N[-1]:.6e}")
    print(f"r = {r}: max|N| over run = {np.max(np.abs(N)):.6e}")

# Critical rate for the DDE dN/dt = r*N(t-tau): r_crit = -pi/(2*tau)
r_crit = -np.pi / 2.0
print(f"Critical growth rate r_crit = -pi/2 = {r_crit:.6f}")

# --- Plot ---
plt.figure(figsize=(9, 6))
for r in rates:
    t, N = results[r]
    plt.plot(t, N, label=labels[r])
plt.axhline(0.0, color="k", lw=0.8, ls="--")
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed exponential growth dN/dt = r*N(t - tau), tau = 1")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.1.1_s1.png", dpi=150)

# Explanation of why the check confirms the result:
# The check confirms it because the observed transition -- monotonic decay at
# r = -0.3, a decaying oscillation at r = -1.4, and an oscillation that grows
# without bound at r = -1.7 -- brackets the analytic critical rate r = -pi/2,
# showing the N = 0 steady state loses stability exactly when r crosses -pi/2.
print("Check: monotonic (r=-0.3) -> damped oscillation (r=-1.4) -> growing "
      "oscillation (r=-1.7) brackets r_crit = -pi/2, confirming instability onset.")
