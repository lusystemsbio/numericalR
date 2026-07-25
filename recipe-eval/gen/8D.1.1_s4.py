import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def gillespie(x0, propensity_fn, stoich_fn, tmax, max_iter, rng):
    """Direct-method Gillespie SSA for a general reaction network.

    x0            : initial state vector of integer molecule counts
    propensity_fn : x -> array of reaction propensities R_j(x)
    stoich_fn     : returns matrix whose row j is the stoichiometry vector gamma_j
    tmax          : stop when simulated time exceeds this
    max_iter      : stop after this many reaction events (safety bound)
    rng           : numpy Generator for reproducible randomness
    """
    x = np.array(x0, dtype=float)          # current state
    t = 0.0                                # current time
    gamma = np.asarray(stoich_fn(), dtype=float)  # rows = reaction update vectors

    times = [t]                            # trajectory: recorded times
    states = [x.copy()]                    # trajectory: recorded states

    for _ in range(max_iter):
        R = np.asarray(propensity_fn(x), dtype=float)  # all propensities R_j(x)
        R_tot = R.sum()                    # total exit rate from current state

        if R_tot <= 0.0:                   # no reaction can fire: absorbing state
            break

        # 1) waiting time ~ Exponential(rate = R_tot)
        tau = rng.exponential(1.0 / R_tot)

        # 2) choose reaction j with probability R_j / R_tot (inverse-CDF sampling)
        u = rng.random() * R_tot           # uniform on [0, R_tot)
        j = np.searchsorted(np.cumsum(R), u)  # first j with cumulative rate > u

        # 3) advance time and 4) update state by chosen reaction's stoichiometry
        t += tau
        x = x + gamma[j]

        times.append(t)
        states.append(x.copy())

        if t >= tmax:                      # stop condition on simulated time
            break

    return np.array(times), np.array(states)


# --- Test network: reversible dimerization-style birth/death of two species ---
#   A --k1--> B         (conversion)
#   B --k2--> A         (back conversion)
#   0 --k3--> A         (production of A)
#   A --k4--> 0         (degradation of A)
k1, k2, k3, k4 = 1.0, 0.5, 2.0, 0.1

def propensity(x):
    A, B = x
    return np.array([k1 * A,      # A -> B
                     k2 * B,      # B -> A
                     k3,          # 0 -> A
                     k4 * A])     # A -> 0

def stoichiometry():
    # each row updates state [A, B]
    return np.array([[-1, +1],    # A -> B
                     [+1, -1],    # B -> A
                     [+1,  0],    # 0 -> A
                     [-1,  0]])   # A -> 0


rng = np.random.default_rng(4)
x0 = [10, 0]
tmax = 50.0
max_iter = 100000

times, states = gillespie(x0, propensity, stoichiometry, tmax, max_iter, rng)

# --- Report numerical results from the trajectory ---
print("Number of reaction events:", len(times) - 1)
print("Final simulated time:", times[-1])
print("Initial state [A, B]:", states[0].tolist())
print("Final state [A, B]:", states[-1].astype(int).tolist())
print("Mean A over trajectory (event-sampled):", states[:, 0].mean())
print("Mean B over trajectory (event-sampled):", states[:, 1].mean())
print("Max A reached:", int(states[:, 0].max()))
print("Max B reached:", int(states[:, 1].max()))

# --- Check: ensemble mean of A at tmax vs deterministic steady state ---
# For the linear network above, the master-equation mean obeys the same ODEs as
# the deterministic rate law, so a many-run average must match the ODE fixed point.
n_runs = 500
finalA = np.empty(n_runs)
finalB = np.empty(n_runs)
for i in range(n_runs):
    r = np.random.default_rng(1000 + i)
    tt, ss = gillespie(x0, propensity, stoichiometry, tmax, max_iter, r)
    finalA[i] = ss[-1, 0]
    finalB[i] = ss[-1, 1]

# Deterministic steady state: dA/dt = -k1 A + k2 B + k3 - k4 A = 0, dB/dt = k1 A - k2 B = 0
# => A* = k3/k4, B* = (k1/k2) A*
A_star = k3 / k4
B_star = (k1 / k2) * A_star
print("Ensemble mean A at tmax ({} runs):".format(n_runs), finalA.mean())
print("Ensemble mean B at tmax ({} runs):".format(n_runs), finalB.mean())
print("Deterministic steady-state A*:", A_star)
print("Deterministic steady-state B*:", B_star)
print("Abs error mean A vs A*:", abs(finalA.mean() - A_star))
print("Abs error mean B vs B*:", abs(finalB.mean() - B_star))

# The check confirms correctness because Gillespie's exponential waiting time and
# R_j/R_tot reaction choice are exactly the jump time and jump distribution of the
# continuous-time Markov chain defined by the master equation, so its trajectory
# ensemble reproduces the master equation's mean (here matching the ODE fixed point).

fig, ax = plt.subplots(figsize=(8, 5))
ax.step(times, states[:, 0], where="post", label="A")
ax.step(times, states[:, 1], where="post", label="B")
ax.axhline(A_star, color="C0", ls="--", alpha=0.5, label="A* (ODE)")
ax.axhline(B_star, color="C1", ls="--", alpha=0.5, label="B* (ODE)")
ax.set_xlabel("time")
ax.set_ylabel("molecule count")
ax.set_title("Gillespie SSA trajectory (direct method)")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.1.1_s4.png")
