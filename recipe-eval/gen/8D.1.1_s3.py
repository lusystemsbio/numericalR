import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def gillespie(state0, propensity_fn, stoich_fn, tmax, max_iter=100000, rng=None):
    """
    Gillespie Stochastic Simulation Algorithm (direct method).

    state0         : initial molecule-count vector (1D array of ints)
    propensity_fn  : function(state) -> 1D array of reaction propensities R_j
    stoich_fn      : function() -> 2D array, row j is the stoichiometry
                     update vector gamma_j applied to the state when
                     reaction j fires
    tmax           : stop when simulated time exceeds this
    max_iter       : hard cap on number of reaction events (stop condition)
    rng            : optional numpy Generator for reproducibility

    Returns (times, states): time points and the state at each point.
    """
    if rng is None:
        rng = np.random.default_rng()

    gamma = np.asarray(stoich_fn(), dtype=float)   # stoichiometry matrix
    state = np.asarray(state0, dtype=float).copy()

    t = 0.0
    times = [t]
    states = [state.copy()]

    for _ in range(max_iter):
        # 1. Compute all propensities R_j for the current state.
        R = np.asarray(propensity_fn(state), dtype=float)
        R_tot = R.sum()

        # If no reaction can fire, the system is frozen -> stop.
        if R_tot <= 0.0:
            break

        # 2. Draw the waiting time tau from an exponential with rate R_tot.
        #    (Time to the next event in a Poisson process of rate R_tot.)
        tau = rng.exponential(1.0 / R_tot)

        # Stop if advancing would pass tmax.
        if t + tau > tmax:
            break

        # 3. Pick which reaction j fires, with probability R_j / R_tot.
        #    Draw u in [0, R_tot) and find the reaction whose cumulative
        #    propensity interval contains u.
        u = rng.random() * R_tot
        j = np.searchsorted(np.cumsum(R), u, side="right")

        # 4. Advance time and update the state by the chosen stoichiometry.
        t += tau
        state = state + gamma[j]

        times.append(t)
        states.append(state.copy())

    return np.asarray(times), np.asarray(states)


# ----------------------------------------------------------------------
# Test network: a simple birth-death process for a single species A.
#   Reaction 1:  0 -> A   (birth)    propensity R1 = k_b            gamma = +1
#   Reaction 2:  A -> 0   (death)    propensity R2 = k_d * n_A      gamma = -1
# This is a generic driver: propensity and stoichiometry are user-supplied.
# ----------------------------------------------------------------------
k_b = 10.0   # birth rate
k_d = 1.0    # per-molecule death rate


def propensity_fn(state):
    n_A = state[0]
    return np.array([k_b, k_d * n_A])


def stoich_fn():
    # Row j is the update vector gamma_j for reaction j.
    return np.array([[+1],    # birth: +1 molecule of A
                     [-1]])   # death: -1 molecule of A


rng = np.random.default_rng(42)
state0 = np.array([0])
tmax = 20.0

times, states = gillespie(state0, propensity_fn, stoich_fn, tmax,
                          max_iter=100000, rng=rng)

n_A = states[:, 0]

print("Number of reaction events (state changes):", len(times) - 1)
print("Final simulated time:", times[-1])
print("Final molecule count n_A:", n_A[-1])
print("Minimum n_A over trajectory:", int(n_A.min()))
print("Maximum n_A over trajectory:", int(n_A.max()))

# ----------------------------------------------------------------------
# Separate check: the SSA generates *exact* samples from the master
# equation, so a long-run / ensemble average of n_A must match the known
# stationary distribution. For the birth-death process the stationary
# distribution is Poisson with mean k_b / k_d, so E[n_A] = k_b / k_d.
#
# One sentence: this check confirms the result because the Gillespie
# direct method draws the exact next-event time and reaction from the
# same waiting-time and selection probabilities defined by the master
# equation, so its sample trajectories are statistically exact and their
# ensemble mean must converge to the master equation's stationary mean.
# ----------------------------------------------------------------------
theory_mean = k_b / k_d

# Ensemble estimate of the stationary mean: run many trajectories and
# average the endpoint state (each endpoint is an independent sample).
n_runs = 2000
endpoints = np.empty(n_runs)
for i in range(n_runs):
    r = np.random.default_rng(1000 + i)
    tt, ss = gillespie(state0, propensity_fn, stoich_fn, tmax,
                       max_iter=100000, rng=r)
    endpoints[i] = ss[-1, 0]

empirical_mean = endpoints.mean()

print("Theoretical stationary mean E[n_A] = k_b/k_d:", theory_mean)
print("Empirical ensemble mean of n_A at t=tmax:", empirical_mean)
print("Absolute difference (empirical - theory):",
      empirical_mean - theory_mean)

# ----------------------------------------------------------------------
# Plot the sample trajectory (step function, since state is piecewise
# constant between events).
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.step(times, n_A, where="post", color="C0", lw=1.2,
        label="SSA trajectory of n_A")
ax.axhline(theory_mean, color="C3", ls="--",
           label="stationary mean k_b/k_d")
ax.set_xlabel("time")
ax.set_ylabel("molecule count n_A")
ax.set_title("Gillespie SSA: birth-death reaction network")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.1.1_s3.png")
print("Saved figure to 8D.1.1_s3.png")
