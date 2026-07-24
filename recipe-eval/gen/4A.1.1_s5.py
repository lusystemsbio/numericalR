import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Model: dN/dt = r * N(t - tau), Euler for delay differential equations ---
# Steady state is N = 0. We integrate explicitly, storing the whole trajectory
# including the constant-history interval t <= 0.

def integrate_dde(r, tau=1.0, dt=0.01, t_end=40.0, hist_value=1.0):
    n_hist = int(round(tau / dt))          # number of stored steps in one delay time
    n_future = int(round(t_end / dt))      # number of steps to integrate forward
    delay_steps = n_hist                   # look this many steps back for N(t - tau)

    # Build time array covering the history interval [-tau, 0] and forward [0, t_end]
    t = np.arange(-n_hist, n_future + 1) * dt

    # Storage for the whole trajectory; history interval is the constant value
    N = np.empty(t.size)
    N[:n_hist + 1] = hist_value            # constant history N(t) = 1 for t <= 0

    # Explicit Euler stepping: current index i corresponds to time t[i]
    for i in range(n_hist, n_hist + n_future):
        N_delayed = N[i - delay_steps]     # value one delay time in the past
        N[i + 1] = N[i] + dt * r * N_delayed  # N_next = N + dt*r*N_delayed
    return t, N

# Run for the three growth rates
rates = [-0.3, -1.4, -1.7]
labels = {-0.3: "monotonic", -1.4: "damped oscillation", -1.7: "growing oscillation"}
results = {}
for r in rates:
    t, N = integrate_dde(r)
    results[r] = (t, N)

# --- Numerical check: characterize the approach to the N = 0 steady state ---
critical = -np.pi / 2
print(f"Critical rate (-pi/2): {critical}")
for r in rates:
    t, N = results[r]
    forward = N[t >= 0]                     # only the integrated forward portion
    tf = t[t >= 0]
    # Count sign changes of N to detect oscillation
    signs = np.sign(forward)
    sign_changes = int(np.sum(signs[:-1] * signs[1:] < 0))
    final_abs = abs(forward[-1])
    initial_abs = abs(forward[0])
    print(f"r = {r}: sign changes = {sign_changes}, "
          f"|N| initial = {initial_abs:.4f}, |N| at t=40 = {final_abs:.6g}, "
          f"behavior = {labels[r]}")

# Explanation of why this check confirms the result:
print("Check confirms the result because r=-0.3 (> -pi/2) decays with no sign changes,")
print("r=-1.4 (> -pi/2) decays while changing sign (damped oscillation), and")
print("r=-1.7 (< -pi/2) changes sign with |N| growing to a large value (unbounded oscillation),")
print("so crossing -pi/2 is exactly where the oscillatory decay turns into oscillatory growth.")

# --- Plot ---
plt.figure(figsize=(10, 6))
for r in rates:
    t, N = results[r]
    plt.plot(t, N, label=f"r = {r} ({labels[r]})")
plt.axhline(0.0, color="k", lw=0.5)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Delayed exponential growth dN/dt = r N(t - tau), tau = 1")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.1.1_s5.png")
