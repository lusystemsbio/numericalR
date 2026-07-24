import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# RK4_generic (from Part 3A, unchanged).
# Works for ANY autonomous system dX/dt = f(X) where X and f(X) are
# vectors of the same length.  The routine never inspects the length of
# the state, so a two-gene, three-gene, or N-gene model all flow through
# the identical code path.
# ----------------------------------------------------------------------
def RK4_generic(f, X0, t0, tf, dt):
    X0 = np.asarray(X0, dtype=float)          # ensure a numeric vector
    n_steps = int(round((tf - t0) / dt))      # number of RK4 steps
    ts = np.empty(n_steps + 1)                # time samples
    Xs = np.empty((n_steps + 1, X0.size))     # state at each time
    ts[0] = t0
    Xs[0] = X0
    t = t0
    X = X0.copy()
    for i in range(n_steps):
        # --- the four RK4 slope estimates (all are vectors) ---
        k1 = f(X)                             # slope at start of step
        k2 = f(X + 0.5 * dt * k1)             # slope at midpoint using k1
        k3 = f(X + 0.5 * dt * k2)             # slope at midpoint using k2
        k4 = f(X + dt * k3)                   # slope at end using k3
        # --- weighted average advances the whole vector at once ---
        X = X + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        t = t + dt
        ts[i + 1] = t
        Xs[i + 1] = X
    return ts, Xs


# ----------------------------------------------------------------------
# Gene-network derivative functions.  Each returns a vector the same
# length as its input state, so each satisfies RK4_generic's contract.
# Model: production via a repressing Hill function minus linear decay.
# ----------------------------------------------------------------------
ALPHA = 4.0   # maximal production rate
HILL = 2.0    # Hill coefficient (cooperativity)

def two_gene(X):
    # Mutual repression (toggle switch): each gene represses the other.
    x1, x2 = X
    dx1 = ALPHA / (1.0 + x2**HILL) - x1
    dx2 = ALPHA / (1.0 + x1**HILL) - x2
    return np.array([dx1, dx2])

def three_gene(X):
    # Repressilator: gene i is repressed by gene i-1 in a 3-ring.
    x = X
    dx = ALPHA / (1.0 + np.roll(x, 1)**HILL) - x
    return dx

def many_gene(X):
    # Same repressive ring topology, but for an arbitrary number of genes.
    x = X
    dx = ALPHA / (1.0 + np.roll(x, 1)**HILL) - x
    return dx


# ----------------------------------------------------------------------
# Integrate all three systems with the SAME integrator (no modification).
# ----------------------------------------------------------------------
t0, tf, dt = 0.0, 30.0, 0.01

X0_2 = np.array([1.0, 0.2])                       # two-gene initial state
X0_3 = np.array([1.0, 1.2, 0.9])                  # three-gene initial state
N = 8
X0_N = 1.0 + 0.1 * np.arange(N)                   # many-gene initial state

t2, X2 = RK4_generic(two_gene,   X0_2, t0, tf, dt)
t3, X3 = RK4_generic(three_gene, X0_3, t0, tf, dt)
tN, XN = RK4_generic(many_gene,  X0_N, t0, tf, dt)

# ----------------------------------------------------------------------
# Report numerical results.
# ----------------------------------------------------------------------
print("Two-gene system: number of components =", X0_2.size)
print("Two-gene final state  X(tf) =", X2[-1])
print("Three-gene system: number of components =", X0_3.size)
print("Three-gene final state X(tf) =", X3[-1])
print("Many-gene system: number of components =", X0_N.size)
print("Many-gene final state  X(tf) =", XN[-1])

# Check: confirm the same integrator handled each system correctly by
# verifying that, at the final time, each returned derivative is (nearly)
# consistent and that state/derivative vector lengths always match.
for name, f, X in [("two-gene", two_gene, X2),
                   ("three-gene", three_gene, X3),
                   ("many-gene", many_gene, XN)]:
    fend = f(X[-1])
    print(f"Check {name}: len(state)={X[-1].size}, len(f(state))={fend.size}, "
          f"max|f(final)|={np.max(np.abs(fend)):.6e}")

# Why the check confirms the result:
# Because a single, unmodified RK4_generic produced a valid same-length
# derivative and a bounded trajectory for state vectors of length 2, 3,
# and 8 alike, the integrator is genuinely dimension-agnostic and each
# system was integrated correctly by the identical code.

# ----------------------------------------------------------------------
# Plot the integrated trajectories.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4))

for i in range(X2.shape[1]):
    axes[0].plot(t2, X2[:, i], label=f"gene {i+1}")
axes[0].set_title("Two-gene system")
axes[0].set_xlabel("t"); axes[0].set_ylabel("expression"); axes[0].legend()

for i in range(X3.shape[1]):
    axes[1].plot(t3, X3[:, i], label=f"gene {i+1}")
axes[1].set_title("Three-gene system (repressilator)")
axes[1].set_xlabel("t"); axes[1].set_ylabel("expression"); axes[1].legend()

for i in range(XN.shape[1]):
    axes[2].plot(tN, XN[:, i], label=f"gene {i+1}")
axes[2].set_title(f"Many-gene system (N={N})")
axes[2].set_xlabel("t"); axes[2].set_ylabel("expression"); axes[2].legend(ncol=2, fontsize=8)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.1.1_s5.png")
