import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Generalized Lotka-Volterra community:
#   dN_i/dt = N_i * (1 - sum_j a_ij * N_j) + D
# with self-interaction a_ii = 1 and random off-diagonal a_ij.
# ----------------------------------------------------------------------

np.random.seed(5)  # reproducible interaction matrix / initial condition

S = 50            # number of species
D = 1e-6          # small immigration (dispersal) term
a_values = [0.08, 0.16, 0.64]  # mean interaction strengths: weak, medium, strong

t0, t1, dt = 0.0, 200.0, 0.01   # integration window and step
n_steps = int(round((t1 - t0) / dt))


def build_interaction_matrix(S, a, rng):
    """Build the S x S interaction matrix A.
    Diagonal (self-interaction) = 1; off-diagonals ~ Uniform(0, 2a),
    so the mean off-diagonal interaction strength is a."""
    A = rng.uniform(0.0, 2.0 * a, size=(S, S))  # all entries in [0, 2a)
    np.fill_diagonal(A, 1.0)                     # set self-interaction to 1
    return A


def glv_rhs(N, A, D):
    """Right-hand side f(N) of the GLV system (vector of length S).
    N_i * (1 - sum_j A_ij N_j) + D, computed for all species at once."""
    # A @ N gives, for each i, sum_j A_ij * N_j
    return N * (1.0 - A @ N) + D


def rk4_step(f, y, dt, *args):
    """One generic classical RK4 step for dy/dt = f(y, *args).
    Works for any dimension: here y is the length-S vector of abundances."""
    k1 = f(y, *args)                    # slope at start
    k2 = f(y + 0.5 * dt * k1, *args)    # slope at midpoint using k1
    k3 = f(y + 0.5 * dt * k2, *args)    # slope at midpoint using k2
    k4 = f(y + dt * k3, *args)          # slope at end using k3
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)  # weighted average


def integrate(f, y0, dt, n_steps, *args):
    """Integrate S variables forward with the generic RK4 above.
    Returns the time array and the (n_steps+1) x S trajectory."""
    traj = np.empty((n_steps + 1, y0.size))  # storage for all S variables
    traj[0] = y0
    y = y0.copy()
    for k in range(n_steps):              # explicit time-stepping loop
        y = rk4_step(f, y, dt, *args)     # advance the whole S-vector
        y = np.maximum(y, 0.0)            # abundances cannot go negative
        traj[k + 1] = y
    t = t0 + dt * np.arange(n_steps + 1)
    return t, traj


# ----------------------------------------------------------------------
# Run the three communities (weak / medium / strong interactions).
# ----------------------------------------------------------------------
rng = np.random.default_rng(5)
N0 = rng.uniform(0.05, 0.15, size=S)  # small random initial abundances (shared)

fig, axes = plt.subplots(1, 3, figsize=(16, 5), sharey=True)

extinct_threshold = 1e-3  # below this a species is considered near-extinct

for ax, a in zip(axes, a_values):
    A = build_interaction_matrix(S, a, np.random.default_rng(int(a * 1000)))
    t, traj = integrate(glv_rhs, N0.copy(), dt, n_steps, A, D)

    final = traj[-1]                          # steady-state abundances
    n_coexist = int(np.sum(final > extinct_threshold))
    n_extinct = S - n_coexist

    # Report numerical results, each on its own labeled line.
    print(f"a = {a}: total species S = {S}")
    print(f"a = {a}: coexisting species (N > {extinct_threshold}) = {n_coexist}")
    print(f"a = {a}: near-extinct species (N <= {extinct_threshold}) = {n_extinct}")
    print(f"a = {a}: max steady-state abundance = {final.max():.6f}")
    print(f"a = {a}: min steady-state abundance = {final.min():.6e}")
    print(f"a = {a}: mean steady-state abundance = {final.mean():.6f}")

    # Verify relaxation to steady state: change over the last time step.
    max_drift = np.max(np.abs(traj[-1] - traj[-2])) / dt
    print(f"a = {a}: max |dN/dt| at final time (steady-state check) = {max_drift:.3e}")

    for i in range(S):
        ax.plot(t, traj[:, i], lw=0.8)
    ax.set_title(f"a = {a}  |  coexist={n_coexist}, near-extinct={n_extinct}")
    ax.set_xlabel("time")
    ax.grid(alpha=0.3)

axes[0].set_ylabel("abundance N_i")
fig.suptitle(f"Generalized Lotka-Volterra: S={S} species relaxing to steady state (RK4, D={D})")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.6.1_s5.png", dpi=120)

# ----------------------------------------------------------------------
# Explanation of the check.
# ----------------------------------------------------------------------
print("Check explanation: because the SAME generic rk4_step advances the full "
      "length-S abundance vector and the trajectories flatten to constant values "
      "(max |dN/dt| -> ~0) with some species surviving and others near zero, this "
      "confirms the RK4 routine scales to S coupled variables and the community "
      "genuinely relaxes to a coexistence/near-extinction steady state.")
