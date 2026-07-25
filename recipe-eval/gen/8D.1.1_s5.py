import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def gillespie(propensity_fn, stoich_fn, x0, tmax, max_iter, rng):
    """
    Gillespie direct method (SSA) for a general reaction network.

    propensity_fn(x) -> array R of reaction rates R_j for current state x
    stoich_fn(j, x)  -> stoichiometry vector gamma_j (state change of reaction j)
    x0               -> initial molecule-count vector
    tmax, max_iter   -> stop conditions
    Returns arrays of times and states (one row per recorded step).
    """
    t = 0.0                       # current time
    x = np.array(x0, dtype=float) # current state (molecule counts)
    times = [t]                   # trajectory storage
    states = [x.copy()]

    for _ in range(max_iter):
        # 1. Compute all propensities and their total rate.
        R = np.asarray(propensity_fn(x), dtype=float)
        R_tot = R.sum()
        if R_tot <= 0.0:          # no reaction can fire -> system frozen
            break

        # 2. Draw the waiting time tau ~ Exponential(rate = R_tot).
        tau = rng.exponential(1.0 / R_tot)

        # 3. Advance time; stop if we would step past tmax.
        if t + tau > tmax:
            break
        t = t + tau

        # 4. Pick reaction j with probability R_j / R_tot (inverse-CDF draw).
        u = rng.random() * R_tot
        j = np.searchsorted(np.cumsum(R), u)

        # 5. Update the state with reaction j's stoichiometry vector.
        x = x + np.asarray(stoich_fn(j, x), dtype=float)

        times.append(t)
        states.append(x.copy())

    return np.array(times), np.array(states)


# ---------------------------------------------------------------------------
# Test network: reversible dimerization-style birth/decay plus conversion.
# Species: [A, B]
# Reactions:
#   0: A -> B          (rate k1 * A)
#   1: B -> A          (rate k2 * B)
#   2: 0 -> A          (rate k3, constant source)
#   3: A -> 0          (rate k4 * A, degradation)
# ---------------------------------------------------------------------------
k1, k2, k3, k4 = 1.0, 0.5, 5.0, 0.2

# Stoichiometry vectors gamma_j (columns = species [A, B]).
STOICH = np.array([
    [-1, +1],   # reaction 0: A -> B
    [+1, -1],   # reaction 1: B -> A
    [+1,  0],   # reaction 2: 0 -> A
    [-1,  0],   # reaction 3: A -> 0
], dtype=float)


def propensity_fn(x):
    A, B = x
    return np.array([k1 * A, k2 * B, k3, k4 * A])


def stoich_fn(j, x):
    return STOICH[j]


rng = np.random.default_rng(5)
x0 = [0.0, 0.0]
tmax = 50.0
max_iter = 100000

times, states = gillespie(propensity_fn, stoich_fn, x0, tmax, max_iter, rng)

# ---------------------------------------------------------------------------
# Reported numerical results from the produced trajectory.
# ---------------------------------------------------------------------------
print("Number of reaction events fired:", len(times) - 1)
print("Final simulated time:", times[-1])
print("Final state A:", states[-1, 0])
print("Final state B:", states[-1, 1])

# Time-average of each species (weighted by the dwell time in each state).
dt = np.diff(times)
if dt.sum() > 0:
    A_timeavg = np.sum(states[:-1, 0] * dt) / dt.sum()
    B_timeavg = np.sum(states[:-1, 1] * dt) / dt.sum()
else:
    A_timeavg = states[0, 0]
    B_timeavg = states[0, 1]
print("Time-averaged A:", A_timeavg)
print("Time-averaged B:", B_timeavg)

# Deterministic steady state for cross-check (mean-field ODE fixed point):
#   dA/dt = -k1 A + k2 B + k3 - k4 A = 0
#   dB/dt =  k1 A - k2 B            = 0  =>  k1 A = k2 B
# From the second equation B = (k1/k2) A; substitute into the first:
#   -k1 A + k1 A + k3 - k4 A = 0  =>  A_ss = k3 / k4
A_ss = k3 / k4
B_ss = (k1 / k2) * A_ss
print("Deterministic steady-state A:", A_ss)
print("Deterministic steady-state B:", B_ss)

# ---------------------------------------------------------------------------
# Independent check: the SSA produces EXACT samples of the master equation.
# For the pure source/degradation subsystem 0 -> A (rate k3), A -> 0 (rate k4*A),
# the master equation has an exact stationary Poisson law with mean lambda = k3/k4.
# We run many independent short trajectories of just that subsystem and compare
# the empirical stationary mean to the exact analytic mean.
# ---------------------------------------------------------------------------
STOICH_bd = np.array([[+1.0], [-1.0]])  # source, degradation on species [A]


def prop_bd(x):
    return np.array([k3, k4 * x[0]])


def stoich_bd(j, x):
    return STOICH_bd[j]


rng2 = np.random.default_rng(12345)
n_trials = 4000
final_A = np.empty(n_trials)
for i in range(n_trials):
    t_bd, s_bd = gillespie(prop_bd, stoich_bd, [0.0], 40.0, 100000, rng2)
    final_A[i] = s_bd[-1, 0]  # sample near stationarity (t=40 >> 1/k4=5)

print("Master-equation check: empirical stationary mean A:", final_A.mean())
print("Master-equation check: empirical stationary var  A:", final_A.var())
print("Master-equation check: exact Poisson mean = var  :", k3 / k4)

# One-sentence explanation:
# Because the SSA draws each inter-event time from the exact exponential
# implied by the total propensity and each reaction with probability
# R_j/R_tot, it samples exactly the same jump process whose probability flow
# IS the chemical master equation, so agreement of the empirical stationary
# distribution with the analytic Poisson law confirms the routine is exact.
print("Note: SSA inter-event times ~ Exp(R_tot) and choice ~ R_j/R_tot are"
      " exactly the master-equation jump process, so matching the analytic"
      " Poisson stationary law confirms exactness.")

# ---------------------------------------------------------------------------
# Plot the time-and-state trajectory (step plot, since states are piecewise
# constant between reaction events).
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5))
ax.step(times, states[:, 0], where="post", label="A")
ax.step(times, states[:, 1], where="post", label="B")
ax.axhline(A_ss, color="C0", ls="--", alpha=0.6, label="A steady state")
ax.axhline(B_ss, color="C1", ls="--", alpha=0.6, label="B steady state")
ax.set_xlabel("time")
ax.set_ylabel("molecule count")
ax.set_title("Gillespie SSA trajectory of the reaction network")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.1.1_s5.png")
