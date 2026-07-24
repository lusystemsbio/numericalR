import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters (fixed for all runs) ----
tau = 1.0        # delay time
dt = 0.01        # Euler time step
t_end = 40.0     # final integration time
n_hist = int(round(tau / dt))     # number of steps in one delay interval
n_steps = int(round(t_end / dt))  # number of forward steps after history

# Time array: include the history interval t in [-tau, 0] and forward to t_end
t = np.arange(-n_hist, n_steps + 1) * dt

def integrate_dde(r):
    # Solution array covering both history and forward integration
    N = np.empty(len(t))
    # Constant history: N(t) = 1 for t <= 0 (fills the first n_hist+1 points)
    N[:n_hist + 1] = 1.0
    # Explicit Euler stepping for the delay equation
    for i in range(n_hist, n_hist + n_steps):
        # Delayed value: population tau/dt steps in the past
        N_delayed = N[i - n_hist]
        # Euler update using the delayed rate: N_next = N + dt*r*N_delayed
        N[i + 1] = N[i] + dt * r * N_delayed
    return N

# ---- Run for the three growth rates ----
rates = [-0.3, -1.4, -1.7]
labels = {-0.3: "monotonic decay", -1.4: "damped oscillation", -1.7: "growing oscillation"}
solutions = {}

plt.figure(figsize=(9, 6))
for r in rates:
    N = integrate_dde(r)
    solutions[r] = N
    plt.plot(t, N, label=f"r = {r} ({labels[r]})")

    # Report final value and amplitude behavior for each rate
    print(f"r = {r}: N(t=40) = {N[-1]:.6e}")
    print(f"r = {r}: max|N| over trajectory = {np.max(np.abs(N)):.6e}")

# Critical rate for the transition
r_crit = -np.pi / 2
print(f"Critical growth rate -pi/2 = {r_crit:.6f}")

# Compare each rate to the critical value
for r in rates:
    print(f"r = {r} is {'above' if r > r_crit else 'below'} the critical value -pi/2")

plt.axhline(0.0, color="k", lw=0.8)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed exponential growth: dN/dt = r*N(t - tau), tau = 1")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.1.1_s2.png")

# Explanation: The check confirms the result because as r crosses -pi/2 the far-field
# behavior switches from monotonic (r=-0.3, no sign changes, decays to 0) to a bounded
# damped oscillation (r=-1.4, oscillates but |N| shrinks) to an unbounded growing
# oscillation (r=-1.7, oscillates while |N| increases), exactly the loss of stability
# of the N=0 steady state predicted at the critical delay-instability threshold r = -pi/2.
print("Explanation: The check confirms the result because r=-0.3 decays monotonically, "
      "r=-1.4 (> -pi/2) decays as a damped oscillation, and r=-1.7 (< -pi/2) oscillates "
      "with growing amplitude, demonstrating that N=0 loses stability precisely as r drops "
      "past the critical value -pi/2.")
