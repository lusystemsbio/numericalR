import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def gillespie(propensity_fn, stoich_matrix, x0, tmax, max_iter=100000, rng=None):
    """
    Gillespie direct-method stochastic simulation algorithm (SSA).

    propensity_fn : function(state) -> 1D array of reaction rates R_j
    stoich_matrix : array of shape (num_reactions, num_species); row j is the
                    stoichiometry vector gamma_j added to the state when
                    reaction j fires
    x0            : initial molecule-count vector
    tmax          : stop when simulation time exceeds this
    max_iter      : stop after at most this many reaction events
    Returns arrays (times, states) giving the exact sample trajectory.
    """
    if rng is None:
        rng = np.random.default_rng()

    stoich_matrix = np.asarray(stoich_matrix, dtype=float)
    x = np.array(x0, dtype=float)
    t = 0.0

    times = [t]
    states = [x.copy()]

    for _ in range(max_iter):
        # 1. compute all propensities R_j at the current state
        R = np.asarray(propensity_fn(x), dtype=float)
        R_tot = R.sum()

        # 2. if no reaction can fire, the system is frozen -> stop
        if R_tot <= 0.0:
            break

        # 3. draw waiting time tau ~ Exponential(rate = R_tot)
        tau = rng.exponential(1.0 / R_tot)
        t = t + tau
        if t > tmax:
            break

        # 4. pick reaction j with probability R_j / R_tot
        j = rng.choice(len(R), p=R / R_tot)

        # 5. advance the state by that reaction's stoichiometry vector
        x = x + stoich_matrix[j]

        times.append(t)
        states.append(x.copy())

    return np.array(times), np.array(states)


# ---------------------------------------------------------------------------
# Test network: reversible dimerization / birth-death style network
#   Species: [A, B]
#   Reaction 0:  A -> B         rate k1 * A
#   Reaction 1:  B -> A         rate k2 * B
#   Reaction 2:  0 -> A         rate k3            (constant production of A)
#   Reaction 3:  A -> 0         rate k4 * A        (degradation of A)
# ---------------------------------------------------------------------------
k1, k2, k3, k4 = 1.0, 0.5, 2.0, 0.3

# stoichiometry: row j = change applied to [A, B] when reaction j fires
gamma = np.array([
    [-1.0, +1.0],   # A -> B
    [+1.0, -1.0],   # B -> A
    [+1.0,  0.0],   # 0 -> A
    [-1.0,  0.0],   # A -> 0
])


def propensities(state):
    A, B = state
    return np.array([k1 * A, k2 * B, k3, k4 * A])


x0 = [10.0, 0.0]
tmax = 20.0
rng = np.random.default_rng(1)

times, states = gillespie(propensities, gamma, x0, tmax, max_iter=100000, rng=rng)

# Report the trajectory summary
print("Number of reaction events:", len(times) - 1)
print("Initial time:", times[0])
print("Final time:", times[-1])
print("Initial state [A, B]:", states[0].tolist())
print("Final state [A, B]:", states[-1].tolist())
print("Mean A over trajectory (event samples):", states[:, 0].mean())
print("Mean B over trajectory (event samples):", states[:, 1].mean())
print("Max A:", states[:, 0].max(), " Min A:", states[:, 0].min())
print("Max B:", states[:, 1].max(), " Min B:", states[:, 1].min())

# ---------------------------------------------------------------------------
# Check: confirm the routine generates EXACT master-equation trajectories.
# For a pure birth-death process  0 -> X (rate a),  X -> 0 (rate b*X),
# the master equation has a known stationary distribution: Poisson(a/b).
# We estimate the stationary mean of X from many independent SSA runs and
# compare it to the exact theoretical mean a/b.
# ---------------------------------------------------------------------------
a, b = 5.0, 1.0            # production and per-molecule degradation rates
gamma_bd = np.array([[+1.0], [-1.0]])   # 0 -> X ; X -> 0


def prop_bd(state):
    X = state[0]
    return np.array([a, b * X])


T_long = 50.0              # run long enough to reach stationarity
n_runs = 400
final_X = np.empty(n_runs)
rng2 = np.random.default_rng(7)
for i in range(n_runs):
    _, s = gillespie(prop_bd, gamma_bd, [0.0], T_long, max_iter=200000, rng=rng2)
    final_X[i] = s[-1, 0]

emp_mean = final_X.mean()
emp_var = final_X.var()
theo_mean = a / b          # Poisson mean
theo_var = a / b           # Poisson variance

print("Birth-death check: empirical stationary mean of X:", emp_mean)
print("Birth-death check: theoretical mean (a/b):", theo_mean)
print("Birth-death check: empirical stationary variance of X:", emp_var)
print("Birth-death check: theoretical variance (a/b):", theo_var)
print("Mean relative error:", abs(emp_mean - theo_mean) / theo_mean)

# Why the check works (one sentence):
print("Explanation: matching the SSA's empirical stationary mean/variance to "
      "the Poisson(a/b) statistics predicted analytically by the master "
      "equation confirms the routine samples exactly from that master equation.")

# ---------------------------------------------------------------------------
# Plot the test-network trajectory
# ---------------------------------------------------------------------------
plt.figure(figsize=(9, 5))
plt.step(times, states[:, 0], where="post", label="A")
plt.step(times, states[:, 1], where="post", label="B")
plt.xlabel("time")
plt.ylabel("molecule count")
plt.title("Gillespie SSA trajectory")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8D.1.1_s1.png")
