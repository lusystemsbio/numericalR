import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def gillespie(x0, propensity_fn, stoich_fn, tmax, max_iter=100000, rng=None):
    """
    Gillespie direct method (stochastic simulation algorithm).

    Parameters
    ----------
    x0            : initial state vector of molecule counts
    propensity_fn : function x -> array of reaction propensities R_j(x)
    stoich_fn     : function j -> stoichiometry vector gamma_j (state update for reaction j)
    tmax          : stop when simulated time exceeds this
    max_iter      : stop after this many reaction events (safety cap)
    rng           : numpy random generator

    Returns
    -------
    times  : 1D array of event times
    states : 2D array (n_events x n_species) of states after each event
    """
    if rng is None:
        rng = np.random.default_rng()

    x = np.array(x0, dtype=float)
    t = 0.0

    # record the trajectory, starting from the initial condition
    times = [t]
    states = [x.copy()]

    for _ in range(max_iter):
        # 1) compute all propensities at the current state
        R = np.asarray(propensity_fn(x), dtype=float)
        R_tot = R.sum()

        # if no reaction can fire, the system is absorbed; stop
        if R_tot <= 0.0:
            break

        # 2) draw the waiting time from an exponential with rate R_tot
        #    tau = -ln(u1)/R_tot  (inverse-transform sample)
        u1 = rng.random()
        tau = -np.log(u1) / R_tot

        # 3) advance time; stop if we would pass tmax
        t = t + tau
        if t > tmax:
            break

        # 4) pick reaction j with probability R_j / R_tot
        #    walk the cumulative sum until it exceeds u2 * R_tot
        u2 = rng.random()
        j = np.searchsorted(np.cumsum(R), u2 * R_tot)

        # 5) update the state by adding reaction j's stoichiometry vector
        x = x + np.asarray(stoich_fn(j), dtype=float)

        times.append(t)
        states.append(x.copy())

    return np.array(times), np.array(states)


# --------------------------------------------------------------------------
# Test network: a simple birth-death / dimerization-style network exercised
# in 8D.2.  Species: [A, B].  Reactions:
#   j=0: A -> B        rate k1 * A          gamma = (-1, +1)
#   j=1: B -> A        rate k2 * B          gamma = (+1, -1)
#   j=2: 0 -> A        rate k3 (birth)      gamma = (+1,  0)
#   j=3: A -> 0        rate k4 * A (death)  gamma = (-1,  0)
# --------------------------------------------------------------------------
k1, k2, k3, k4 = 1.0, 0.5, 2.0, 0.3

STOICH = np.array([[-1, +1],
                   [+1, -1],
                   [+1,  0],
                   [-1,  0]], dtype=float)


def propensities(x):
    A, B = x
    return np.array([k1 * A, k2 * B, k3, k4 * A])


def stoichiometry(j):
    return STOICH[j]


rng = np.random.default_rng(42)
x0 = [10.0, 0.0]
tmax = 20.0

times, states = gillespie(x0, propensities, stoichiometry,
                          tmax=tmax, max_iter=100000, rng=rng)

# Report trajectory summary
print("Number of reaction events:", len(times) - 1)
print("Final time:", times[-1])
print("Final state [A, B]:", states[-1])
print("Initial state [A, B]:", states[0])
print("Max A over trajectory:", states[:, 0].max())
print("Max B over trajectory:", states[:, 1].max())
print("Time-averaged A:",
      np.trapz(states[:-1, 0], times[:-1]) / times[-2] if len(times) > 1 else states[0, 0])
print("Time-averaged B:",
      np.trapz(states[:-1, 1], times[:-1]) / times[-2] if len(times) > 1 else states[0, 1])

# --------------------------------------------------------------------------
# Check: confirm the routine generates EXACT master-equation trajectories.
# Take the pure birth-death of a single species A with birth rate a and
# death rate b*A.  Its stationary distribution is Poisson(a/b).  Averaging
# many long Gillespie runs must reproduce that stationary mean a/b and the
# Poisson property mean == variance.
# --------------------------------------------------------------------------
a, b = 5.0, 1.0
STOICH_BD = np.array([[+1], [-1]], dtype=float)


def prop_bd(x):
    return np.array([a, b * x[0]])


def stoich_bd(j):
    return STOICH_BD[j]


rng2 = np.random.default_rng(7)
n_samples = 4000
sample_time = 50.0  # long enough to reach stationarity
samples = []
for _ in range(n_samples):
    tt, ss = gillespie([0.0], prop_bd, stoich_bd,
                       tmax=sample_time, max_iter=100000, rng=rng2)
    samples.append(ss[-1, 0])  # state sampled at the final time
samples = np.array(samples)

print("Birth-death check: theoretical stationary mean (a/b):", a / b)
print("Birth-death check: empirical mean:", samples.mean())
print("Birth-death check: empirical variance (Poisson => equals mean):", samples.var())

# One-sentence explanation of why this check works:
print("Why the check confirms the result: because the direct method draws "
      "the exact exponential waiting time and exact reaction probabilities of "
      "the master equation, its samples must converge to the master equation's "
      "exact stationary distribution (here Poisson(a/b) with mean == variance).")

# --------------------------------------------------------------------------
# Plot the two-species test trajectory as a step function.
# --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 5))
ax.step(times, states[:, 0], where="post", label="A")
ax.step(times, states[:, 1], where="post", label="B")
ax.set_xlabel("time")
ax.set_ylabel("molecule count")
ax.set_title("Gillespie direct-method trajectory")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.1.1_s2.png")
