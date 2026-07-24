import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model: dN/dt = r * N(t - tau), steady state N = 0 ----
# Euler for a delay differential equation:
#   N_next = N + dt * r * N_delayed
# where N_delayed is the value tau/dt integration steps in the past.

# ---- Parameters ----
tau = 1.0            # delay time
dt = 0.01            # Euler time step
t_end = 40.0         # final integration time
r_values = [-0.3, -1.4, -1.7]   # growth rates to test

delay_steps = int(round(tau / dt))   # how many steps back the delayed value sits (=100)

# ---- Time grid: include the history interval [-tau, 0] then [0, t_end] ----
# Index 0 corresponds to t = -tau; the delayed lookup is always index i - delay_steps.
t = np.arange(-tau, t_end + dt, dt)
n_total = len(t)
zero_index = delay_steps   # index where t = 0 (end of history)

print(f"delay_steps (tau/dt) = {delay_steps}")
print(f"number of stored points (incl. history) = {n_total}")
print(f"index of t=0 in trajectory = {zero_index}")

results = {}

for r in r_values:
    # Allocate the whole trajectory, history interval included.
    N = np.zeros(n_total)

    # Constant history: N(t) = 1 for t <= 0  (fills the history interval and t = 0).
    N[:zero_index + 1] = 1.0

    # Explicit Euler stepping forward from t = 0.
    for i in range(zero_index, n_total - 1):
        N_delayed = N[i - delay_steps]     # value one delay time in the past
        N[i + 1] = N[i] + dt * r * N_delayed   # Euler update

    results[r] = N

    # Report a few diagnostic numbers per r.
    N_final = N[-1]
    N_min = N.min()
    N_max = N.max()
    print(f"r = {r}: N(t=40) = {N_final:.6g}, min(N) = {N_min:.6g}, max(N) = {N_max:.6g}")

# ---- Plot ----
plt.figure(figsize=(9, 6))
for r in r_values:
    plt.plot(t, results[r], label=f"r = {r}")
plt.axhline(0.0, color="gray", lw=0.8, ls="--")
plt.axvline(0.0, color="gray", lw=0.5)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed exponential growth: dN/dt = r N(t - tau), tau = 1")
plt.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.1.1_s4.png")

# ---- Stability check across the critical value r_crit = -pi/2 ----
r_crit = -np.pi / 2
print(f"critical rate r_crit = -pi/2 = {r_crit:.6f}")

# Classify each trajectory by counting sign changes (oscillation) and
# by whether the late-time envelope grows or decays.
for r in r_values:
    N = results[r]
    post = N[zero_index:]                      # part from t = 0 onward
    sign_changes = int(np.sum(np.diff(np.sign(post)) != 0))
    early_amp = np.max(np.abs(post[:len(post)//2]))
    late_amp = np.max(np.abs(post[len(post)//2:]))
    growing = late_amp > early_amp
    if sign_changes == 0:
        behavior = "monotonic decay"
    elif growing:
        behavior = "growing oscillation (unbounded)"
    else:
        behavior = "damped oscillation"
    print(f"r = {r}: sign changes = {sign_changes}, "
          f"early amp = {early_amp:.4g}, late amp = {late_amp:.4g} -> {behavior}")

# ---- One-sentence explanation ----
print("Explanation: As r decreases past -pi/2, the dominant root of the "
      "characteristic equation lambda = r*exp(-lambda*tau) moves from real-negative "
      "(monotonic decay) to complex with negative real part (damped oscillation) to "
      "complex with positive real part (growing oscillation), so observing decay change "
      "from monotonic (r=-0.3) to damped (r=-1.4) to unbounded growth (r=-1.7) confirms "
      "the critical crossing of the N=0 steady state's stability at r=-pi/2.")
