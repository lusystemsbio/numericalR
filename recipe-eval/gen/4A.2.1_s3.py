import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Problem setup: delayed exponential growth  dN/dt = r * N(t - tau) ----
r   = -1.7          # unstable case
tau = 1.0           # delay
dt  = 0.01          # time step
T   = 40.0          # final time
d   = int(round(tau / dt))   # delay expressed in whole time steps (= 100)

n_steps = int(round(T / dt))
t = np.linspace(0.0, T, n_steps + 1)

# Constant history N(t) = 1 for t <= 0.
# We store the solution in an array long enough to also index the history:
# index i corresponds to time (i - d) * dt, so indices 0..d-1 hold the history.
def rhs(Ndelayed):
    # right-hand side depends only on the delayed value
    return r * Ndelayed

# ---------- First-order Euler ----------
Ne = np.ones(n_steps + 1 + d)   # first d entries = history value 1
for i in range(d, d + n_steps):
    # delayed value N(t - tau) is the stored value one delay-length back
    N_delay = Ne[i - d]
    Ne[i + 1] = Ne[i] + dt * rhs(N_delay)
euler = Ne[d:]   # drop the history padding

# ---------- Second-order Heun (Euler predictor + trapezoidal corrector) ----------
Nh = np.ones(n_steps + 1 + d)
for i in range(d, d + n_steps):
    N_delay      = Nh[i - d]        # delayed value at current time  t_i
    N_delay_next = Nh[i + 1 - d]    # delayed value at next time     t_{i+1}
    #                                 (already stored, since delay = whole steps)

    # 1) Euler predictor: provisional step forward
    N_pred = Nh[i] + dt * rhs(N_delay)

    # 2) Trapezoidal corrector: average the derivative at t_i and t_{i+1}.
    #    Both derivatives use only stored delayed values, so the predictor
    #    supplies the forward point while the corrector refines the slope.
    Nh[i + 1] = Nh[i] + 0.5 * dt * (rhs(N_delay) + rhs(N_delay_next))
heun = Nh[d:]   # drop the history padding

# ---------- Numerical results ----------
print(f"Delay in whole steps d = {d}")
print(f"Euler final value N({T}) = {euler[-1]:.6f}")
print(f"Heun  final value N({T}) = {heun[-1]:.6f}")

# Early-vs-late comparison
i_early = int(round(2.0 / dt))    # t = 2
i_late  = int(round(30.0 / dt))   # t = 30
print(f"Euler at t=2  = {euler[i_early]:.6f}")
print(f"Heun  at t=2  = {heun[i_early]:.6f}")
print(f"|Euler-Heun| at t=2  = {abs(euler[i_early]-heun[i_early]):.6e}")
print(f"Euler at t=30 = {euler[i_late]:.6f}")
print(f"Heun  at t=30 = {heun[i_late]:.6f}")
print(f"|Euler-Heun| at t=30 = {abs(euler[i_late]-heun[i_late]):.6e}")

max_abs_diff = np.max(np.abs(euler - heun))
print(f"Max |Euler - Heun| over [0,{T}] = {max_abs_diff:.6f}")

# ---------- Consistency check ----------
early_agree = abs(euler[i_early] - heun[i_early]) < 1e-2
drift_apart = abs(euler[i_late] - heun[i_late]) > 10.0 * abs(euler[i_early] - heun[i_early] + 1e-30)
euler_overshoots = np.max(np.abs(euler)) > np.max(np.abs(heun))
print(f"Check - trajectories agree early (t=2): {early_agree}")
print(f"Check - trajectories drift apart later: {drift_apart}")
print(f"Check - Euler curve overshoots Heun (larger amplitude): {euler_overshoots}")
print(f"Euler peak |N| = {np.max(np.abs(euler)):.6f}")
print(f"Heun  peak |N| = {np.max(np.abs(heun)):.6f}")

# The check confirms the result because a first-order method accumulates more
# truncation error per step than a second-order one, so identical early behavior
# followed by growing separation with Euler's larger swings is exactly the
# signature of Euler being the less accurate (over-amplifying) scheme.

# ---------- Plot ----------
plt.figure(figsize=(10, 6))
plt.plot(t, euler, label="Euler (1st order)", lw=1.5)
plt.plot(t, heun,  label="Heun (2nd order)", lw=1.5, ls="--")
plt.axhline(0.0, color="k", lw=0.5)
plt.xlabel("t")
plt.ylabel("N(t)")
plt.title(f"Delayed growth dN/dt = r N(t-tau), r={r}, tau={tau}, dt={dt}")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4A.2.1_s3.png")
print("Saved plot to 4A.2.1_s3.png")
