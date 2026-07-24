import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.1.1_s3.png"

# ----------------------------------------------------------------------
# Generic RK4 from Part 3A, UNCHANGED.
# It works for any autonomous system dX/dt = f(X) as long as f returns a
# vector the same length as the state vector X.  It never looks at the
# dimension of X, so the same routine serves every system below.
# ----------------------------------------------------------------------
def RK4_generic(f, X0, t):
    X0 = np.asarray(X0, dtype=float)          # state vector (any length)
    X = np.zeros((len(t), len(X0)))           # storage: one row per time
    X[0] = X0
    for i in range(len(t) - 1):
        h = t[i + 1] - t[i]                   # step size
        Xi = X[i]
        k1 = f(Xi)                            # slope at start
        k2 = f(Xi + 0.5 * h * k1)             # slope at midpoint (using k1)
        k3 = f(Xi + 0.5 * h * k2)             # slope at midpoint (using k2)
        k4 = f(Xi + h * k3)                   # slope at end
        X[i + 1] = Xi + (h / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)  # weighted avg
    return X

# ----------------------------------------------------------------------
# Gene-regulatory models.  Each returns a vector f(X) the same length as X.
# Repression uses a Hill function; each gene represses the "next" one in a
# ring, so the identical mathematical form scales to any number of genes.
# ----------------------------------------------------------------------
def hill_repression(x, K=1.0, n=2.0):
    return 1.0 / (1.0 + (x / K) ** n)         # high repressor -> low output

# Two-gene mutual-repression toggle switch: gene0 <-> gene1
def f_two_gene(X, alpha=4.0, beta=1.0):
    x0, x1 = X
    return np.array([
        alpha * hill_repression(x1) - beta * x0,   # gene0 repressed by gene1
        alpha * hill_repression(x0) - beta * x1,   # gene1 repressed by gene0
    ])

# Three-gene repressilator ring: 0 -| 1 -| 2 -| 0
def f_three_gene(X, alpha=4.0, beta=1.0):
    x = np.asarray(X)
    N = len(x)
    prev = x[(np.arange(N) - 1) % N]           # repressor of each gene
    return alpha * hill_repression(prev) - beta * x

# Many-gene ring repressilator: gene i repressed by gene i-1 (mod N)
def f_many_gene(X, alpha=4.0, beta=1.0):
    x = np.asarray(X)
    N = len(x)
    prev = x[(np.arange(N) - 1) % N]           # ring topology, any N
    return alpha * hill_repression(prev) - beta * x

# ----------------------------------------------------------------------
# Integrate each system with the SAME RK4_generic (no modification).
# ----------------------------------------------------------------------
t = np.linspace(0.0, 40.0, 4001)              # shared time grid

X2 = RK4_generic(f_two_gene,   [1.0, 0.2],                 t)   # 2 genes
X3 = RK4_generic(f_three_gene, [1.0, 0.5, 0.2],            t)   # 3 genes
N_many = 8
X0_many = 0.5 + 0.5 * np.cos(np.arange(N_many))            # spread initial conditions
Xm = RK4_generic(f_many_gene,  X0_many,                    t)   # N genes

# ----------------------------------------------------------------------
# Print numerical results (final states and a couple of diagnostics).
# ----------------------------------------------------------------------
print("Two-gene system: state length =", X2.shape[1])
for i in range(X2.shape[1]):
    print(f"  gene{i} final value = {X2[-1, i]:.6f}")

print("Three-gene system: state length =", X3.shape[1])
for i in range(X3.shape[1]):
    print(f"  gene{i} final value = {X3[-1, i]:.6f}")

print("Many-gene system: state length =", Xm.shape[1])
for i in range(Xm.shape[1]):
    print(f"  gene{i} final value = {Xm[-1, i]:.6f}")

# The check: the integrator handled lengths 2, 3, and N with one routine.
lengths_handled = [X2.shape[1], X3.shape[1], Xm.shape[1]]
print("State lengths handled by the single RK4_generic:", lengths_handled)
print("All finite (no blow-up)?", bool(np.all(np.isfinite(X2)) and
                                       np.all(np.isfinite(X3)) and
                                       np.all(np.isfinite(Xm))))
print("Three-gene oscillation amplitude (gene0, last half) =",
      f"{X3[len(t)//2:, 0].max() - X3[len(t)//2:, 0].min():.6f}")
print("Many-gene oscillation amplitude (gene0, last half) =",
      f"{Xm[len(t)//2:, 0].max() - Xm[len(t)//2:, 0].min():.6f}")

# ----------------------------------------------------------------------
# Plot the three trajectories.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True)

for i in range(X2.shape[1]):
    axes[0].plot(t, X2[:, i], label=f"gene{i}")
axes[0].set_title("Two-gene toggle switch")
axes[0].set_ylabel("expression")
axes[0].legend(loc="upper right")

for i in range(X3.shape[1]):
    axes[1].plot(t, X3[:, i], label=f"gene{i}")
axes[1].set_title("Three-gene repressilator")
axes[1].set_ylabel("expression")
axes[1].legend(loc="upper right")

for i in range(Xm.shape[1]):
    axes[2].plot(t, Xm[:, i], label=f"gene{i}")
axes[2].set_title(f"Many-gene ({N_many}) ring repressilator")
axes[2].set_ylabel("expression")
axes[2].set_xlabel("time")
axes[2].legend(loc="upper right", ncol=2, fontsize=8)

fig.tight_layout()
plt.savefig(OUT)

# One-sentence explanation of why the check confirms the result:
print("Explanation: Because the identical, unmodified RK4_generic produced "
      "finite, sensible trajectories for state vectors of length 2, 3, and "
      f"{N_many}, it confirms the integrator is truly dimension-agnostic and "
      "depends only on f(X) returning a vector matching the state length.")
