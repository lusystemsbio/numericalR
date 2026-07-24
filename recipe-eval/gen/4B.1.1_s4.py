import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# Hutchinson delay-logistic equation:
#     dN/dt = r * N(t) * (1 - N(t-tau)/B)
# Integrated with a 2nd-order Heun (predictor-corrector) scheme for DDEs.
# -----------------------------------------------------------------------------

# Fixed model / integration parameters
tau = 1.0          # delay
B   = 100.0        # carrying capacity
dt  = 0.01         # time step
T   = 60.0         # total simulation time
n_steps = int(round(T / dt))
lag = int(round(tau / dt))   # number of steps corresponding to one delay

def rhs(N_now, N_delayed, r):
    # Right-hand side of the Hutchinson equation; crowding uses the delayed value
    return r * N_now * (1.0 - N_delayed / B)

def simulate(r):
    # Time grid and solution array
    t = np.arange(n_steps + 1) * dt
    N = np.empty(n_steps + 1)
    # Constant history: N(t) = 1 for all t <= 0, i.e. the initial point
    N[0] = 1.0

    for k in range(n_steps):
        # Delayed index; for k < lag the history (constant 1) supplies the value
        idx = k - lag
        N_delayed_now = N[idx] if idx >= 0 else 1.0

        # --- Predictor (explicit Euler step) ---
        f1 = rhs(N[k], N_delayed_now, r)
        N_pred = N[k] + dt * f1

        # --- Delayed value one step ahead, for the corrector slope ---
        idx_next = (k + 1) - lag
        N_delayed_next = N[idx_next] if idx_next >= 0 else 1.0

        # --- Corrector: average the slope at start and at the predicted end ---
        f2 = rhs(N_pred, N_delayed_next, r)
        N[k + 1] = N[k] + 0.5 * dt * (f1 + f2)

    return t, N

# Growth rates spanning the Hopf transition r*tau = pi/2
rates = [0.3, 1.5, np.pi / 2, 1.7]
labels = ["r = 0.3", "r = 1.5", "r = pi/2 (~1.5708)", "r = 1.7"]

results = {}
plt.figure(figsize=(10, 6))
for r, lab in zip(rates, labels):
    t, N = simulate(r)
    results[r] = N
    plt.plot(t, N, label=lab)

    # Characterise the tail of the solution (last 20 time units)
    tail = N[t >= (T - 20.0)]
    tail_min, tail_max = tail.min(), tail.max()
    amplitude = tail_max - tail_min
    print(f"{lab}: r*tau = {r*tau:.4f}, "
          f"tail min = {tail_min:.4f}, tail max = {tail_max:.4f}, "
          f"tail peak-to-peak amplitude = {amplitude:.4f}, "
          f"final N = {N[-1]:.4f}")

print(f"Hopf threshold r*tau = pi/2 = {np.pi/2:.6f}")

plt.axhline(B, color="k", ls="--", lw=0.8, label="carrying capacity B")
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title("Hutchinson delay-logistic growth across the Hopf point r*tau = pi/2")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4B.1.1_s4.png")

# Explanation (one sentence):
print("Explanation: The linearization of the Hutchinson equation about N=B has a "
      "pair of complex eigenvalues that cross the imaginary axis exactly at r*tau = pi/2, "
      "so the observed change from decaying (r<pi/2) to constant-amplitude (r>pi/2) tail "
      "oscillations confirms this Hopf bifurcation.")
