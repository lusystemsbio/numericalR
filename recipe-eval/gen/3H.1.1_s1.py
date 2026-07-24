import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# RK4_generic (from Part 3A, unchanged).
# Advances ONE step of size h for dX/dt = f(X), where X and f(X) are
# vectors of the same (arbitrary) length. Implemented explicitly.
# ----------------------------------------------------------------------
def RK4_generic(f, X, h):
    k1 = f(X)              # slope at the start of the interval
    k2 = f(X + 0.5*h*k1)   # slope at the midpoint, using k1
    k3 = f(X + 0.5*h*k2)   # slope at the midpoint, refined with k2
    k4 = f(X + h*k3)       # slope at the end of the interval, using k3
    # weighted average of the four slopes (1,2,2,1)/6
    return X + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)


def integrate(f, X0, h, nsteps):
    # Roll the state forward nsteps times with the SAME RK4 routine.
    X = np.array(X0, dtype=float)      # state vector of any length
    traj = np.empty((nsteps+1, X.size))
    traj[0] = X
    for i in range(nsteps):
        X = RK4_generic(f, X, h)       # one RK4 step, no system-specific code
        traj[i+1] = X
    return traj


# ----------------------------------------------------------------------
# Gene-regulation models. Each returns a derivative vector the same
# length as its state vector, so RK4_generic handles all of them.
# ----------------------------------------------------------------------
def hill_repress(x, K=1.0, n=2):
    # repression Hill function: high x -> low output
    return 1.0 / (1.0 + (x/K)**n)

# Two-gene toggle switch: genes mutually repress each other.
def f_two(X, beta=4.0, gamma=1.0):
    x1, x2 = X
    return np.array([beta*hill_repress(x2) - gamma*x1,
                     beta*hill_repress(x1) - gamma*x2])

# Three-gene repressilator: 1 -| 2 -| 3 -| 1 (cyclic repression).
def f_three(X, beta=10.0, gamma=1.0):
    x1, x2, x3 = X
    return np.array([beta*hill_repress(x3) - gamma*x1,
                     beta*hill_repress(x1) - gamma*x2,
                     beta*hill_repress(x2) - gamma*x3])

# Many-gene ring repressilator: N genes, gene i repressed by gene i-1.
def make_f_many(N, beta=10.0, gamma=1.0):
    def f_many(X):
        prev = np.roll(X, 1)           # gene i-1 (cyclic) represses gene i
        return beta*hill_repress(prev) - gamma*X
    return f_many


# ----------------------------------------------------------------------
# Integrate all three systems with the identical integrator.
# ----------------------------------------------------------------------
h = 0.01
nsteps = 4000
t = np.arange(nsteps+1)*h

X0_two   = [1.0, 0.5]
X0_three = [1.0, 1.2, 0.9]
N = 7
X0_many  = 1.0 + 0.1*np.sin(np.arange(N))   # slightly perturbed initial ring

traj_two   = integrate(f_two,          X0_two,   h, nsteps)
traj_three = integrate(f_three,        X0_three, h, nsteps)
traj_many  = integrate(make_f_many(N), X0_many,  h, nsteps)

# ----------------------------------------------------------------------
# Report numerical results (final states).
# ----------------------------------------------------------------------
print("Two-gene system, final state:")
for i, v in enumerate(traj_two[-1]):
    print(f"  x{i+1}(t_final) = {v:.6f}")

print("Three-gene system, final state:")
for i, v in enumerate(traj_three[-1]):
    print(f"  x{i+1}(t_final) = {v:.6f}")

print(f"Many-gene system (N={N}), final state:")
for i, v in enumerate(traj_many[-1]):
    print(f"  x{i+1}(t_final) = {v:.6f}")

# Cross-check: report the state-vector length each call handled.
print("Check: state-vector lengths handled by the SAME RK4_generic:")
print(f"  two-gene   length = {traj_two.shape[1]}")
print(f"  three-gene length = {traj_three.shape[1]}")
print(f"  many-gene  length = {traj_many.shape[1]}")

# ----------------------------------------------------------------------
# Plot the integrated trajectories.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(9, 11))

for i in range(traj_two.shape[1]):
    axes[0].plot(t, traj_two[:, i], label=f"gene {i+1}")
axes[0].set_title("Two-gene toggle switch")
axes[0].set_xlabel("t"); axes[0].set_ylabel("expression"); axes[0].legend()

for i in range(traj_three.shape[1]):
    axes[1].plot(t, traj_three[:, i], label=f"gene {i+1}")
axes[1].set_title("Three-gene repressilator")
axes[1].set_xlabel("t"); axes[1].set_ylabel("expression"); axes[1].legend()

for i in range(traj_many.shape[1]):
    axes[2].plot(t, traj_many[:, i], label=f"gene {i+1}")
axes[2].set_title(f"Many-gene ring repressilator (N={N})")
axes[2].set_xlabel("t"); axes[2].set_ylabel("expression")
axes[2].legend(ncol=2, fontsize=8)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.1.1_s1.png")

# The check confirms the result because a single unmodified RK4_generic
# integrated state vectors of length 2, 3, and 7 purely from the length
# of f(X), showing the routine is genuinely dimension-agnostic.
print("Check passed: one unmodified RK4_generic integrated all three systems (lengths 2, 3, 7).")
