import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Hutchinson delay-logistic (delay-logistic) equation:
#   dN/dt = r * N(t) * (1 - N(t - tau) / B)
# The crowding term uses the population one delay tau in the past.
# Integrated with a 2nd-order Heun (predictor-corrector) DDE scheme.
# ---------------------------------------------------------------

# Model / simulation parameters
tau = 1.0          # delay
B   = 100.0        # carrying capacity
dt  = 0.01         # time step
Tmax = 200.0       # total simulated time (long enough to see limit cycle)
N0  = 1.0          # constant history value:  N(t) = 1 for t <= 0

d = int(round(tau / dt))          # delay expressed in integer time steps (=100)
nsteps = int(round(Tmax / dt))    # number of integration steps
t = np.linspace(0.0, Tmax, nsteps + 1)

# right-hand side f = dN/dt, depends on current N and the delayed N
def f(N, N_delayed, r):
    return r * N * (1.0 - N_delayed / B)

def integrate(r):
    N = np.empty(nsteps + 1)
    N[0] = N0
    for j in range(nsteps):
        # delayed value at current time t_j  (history = N0 before t=0)
        N_del_now  = N[j - d] if (j - d) >= 0 else N0
        # delayed value at next time t_{j+1} (needed by the corrector)
        N_del_next = N[j + 1 - d] if (j + 1 - d) >= 0 else N0

        # --- Heun step ---
        k1 = f(N[j], N_del_now, r)              # slope at start
        N_pred = N[j] + dt * k1                 # Euler predictor
        k2 = f(N_pred, N_del_next, r)           # slope at predicted end point
        N[j + 1] = N[j] + 0.5 * dt * (k1 + k2)  # averaged (corrected) step
    return N

# Growth rates spanning the Hopf threshold r*tau = pi/2
rs = [0.3, 1.5, np.pi / 2.0, 1.7]
labels = ["r = 0.3", "r = 1.5", "r = pi/2 (~1.5708)", "r = 1.7"]
sols = {}

print(f"Hopf threshold value r*tau = pi/2 = {np.pi/2:.6f}")
print(f"delay steps d = tau/dt = {d}")
print("")

for r, lab in zip(rs, labels):
    N = integrate(r)
    sols[r] = N
    # analyse the tail (last 40 time units) to classify behaviour
    tail = N[t >= (Tmax - 40.0)]
    amp = 0.5 * (tail.max() - tail.min())   # half peak-to-peak = oscillation amplitude
    print(f"{lab}:  r*tau = {r*tau:.4f}")
    print(f"    final N(Tmax)         = {N[-1]:.6f}")
    print(f"    tail mean             = {tail.mean():.6f}")
    print(f"    tail oscillation amp  = {amp:.6f}")
    print("")

# ---------------------------------------------------------------
# Plot N(t) for each r
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))
for r, lab in zip(rs, labels):
    ax.plot(t, sols[r], label=lab)
ax.axhline(B, color="k", ls="--", lw=0.8, label="carrying capacity B")
ax.set_xlabel("t")
ax.set_ylabel("N(t)")
ax.set_title("Hutchinson delay-logistic growth (Heun DDE integrator)\n"
             "from monotonic growth to a sustained limit cycle")
ax.legend(loc="best")
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/4B.1.1_s3.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the check confirms the result:
# ---------------------------------------------------------------
print("Explanation: Linear stability of the equilibrium N=B loses stability exactly")
print("when r*tau = pi/2, so r=0.3 and r=1.5 (r*tau < pi/2) return to B (monotonic,")
print("then damped-oscillatory), r=pi/2 sits at the neutral Hopf point and decays")
print("only marginally, and r=1.7 (r*tau > pi/2) grows into a stable limit cycle —")
print("this ordered transition across r*tau = pi/2 confirms the Hopf bifurcation.")
