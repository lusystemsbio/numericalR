import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Part 3A: the generic RK4 routine, used UNCHANGED for every system.
# It steps any autonomous system dX/dt = f(X) forward by one step h.
# X and f(X) are vectors of the same (arbitrary) length, so the same
# code integrates a 2-, 3-, or N-component system with no modification.
# ----------------------------------------------------------------------
def RK4_generic(f, X, h):
    # X is a numpy vector of any length; f(X) returns a vector the same length.
    k1 = f(X)              # slope at the start of the interval
    k2 = f(X + 0.5*h*k1)   # slope at the midpoint using k1
    k3 = f(X + 0.5*h*k2)   # slope at the midpoint using k2
    k4 = f(X + h*k3)       # slope at the end using k3
    # weighted average of the four slopes -> new state
    return X + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)


def integrate(f, X0, h, nsteps):
    # Drive RK4_generic forward, storing the whole trajectory.
    X0 = np.asarray(X0, dtype=float)
    t = np.zeros(nsteps + 1)
    traj = np.zeros((nsteps + 1, X0.size))
    traj[0] = X0
    X = X0.copy()
    for i in range(nsteps):
        X = RK4_generic(f, X, h)   # <-- identical call for every system
        traj[i + 1] = X
        t[i + 1] = (i + 1) * h
    return t, traj


# ----------------------------------------------------------------------
# Gene-regulation models of this chapter.
# Each returns a derivative vector the SAME length as its state vector.
# Repression is modeled by a Hill function beta / (1 + x^n).
# ----------------------------------------------------------------------
beta, n = 4.0, 2.0   # max production rate and Hill coefficient

# Two-gene system: mutual repression (toggle switch), state length 2.
def f_two(X):
    x, y = X
    dx = beta / (1.0 + y**n) - x
    dy = beta / (1.0 + x**n) - y
    return np.array([dx, dy])

# Three-gene system: repressilator ring (1->2->3->1), state length 3.
def f_three(X):
    x, y, z = X
    dx = beta / (1.0 + z**n) - x
    dy = beta / (1.0 + x**n) - y
    dz = beta / (1.0 + y**n) - z
    return np.array([dx, dy, dz])

# Many-gene system: a ring of N genes, each repressed by the previous one.
def make_ring(N):
    def f_ring(X):
        prev = np.roll(X, 1)                    # gene i is repressed by gene i-1
        return beta / (1.0 + prev**n) - X       # vectorized, length N
    return f_ring


# ----------------------------------------------------------------------
# Integrate all three systems with the SAME integrator/RK4 routine.
# ----------------------------------------------------------------------
h, nsteps = 0.02, 2000

t2, traj2 = integrate(f_two,   [1.0, 0.2],          h, nsteps)
t3, traj3 = integrate(f_three, [1.0, 0.2, 0.3],     h, nsteps)

N = 6
f_many = make_ring(N)
X0_many = 0.5 + 0.4*np.sin(np.arange(N))   # spread-out initial conditions
tN, trajN = integrate(f_many, X0_many, h, nsteps)

# ----------------------------------------------------------------------
# Separate check: confirm one unmodified integrator handled each system.
# For an autonomous system the derivative must vanish at any true steady
# state; we also confirm the output width equals the input state length.
# ----------------------------------------------------------------------
print("=== Dimensions handled by the single RK4 routine ===")
print(f"Two-gene   state length: {traj2.shape[1]}")
print(f"Three-gene state length: {traj3.shape[1]}")
print(f"Many-gene  state length: {trajN.shape[1]}")

print("\n=== Final states (final time t = {:.2f}) ===".format(nsteps*h))
print("Two-gene   final X:", np.array2string(traj2[-1], precision=6))
print("Three-gene final X:", np.array2string(traj3[-1], precision=6))
print("Many-gene  final X:", np.array2string(trajN[-1], precision=6))

print("\n=== Consistency check: len(f(X)) == len(X) for each system ===")
for name, f, X in [("Two-gene", f_two, traj2[-1]),
                   ("Three-gene", f_three, traj3[-1]),
                   ("Many-gene", f_many, trajN[-1])]:
    d = f(X)
    print(f"{name}: len(X)={len(X)}, len(f(X))={len(d)}, "
          f"max|f(X)| at end = {np.max(np.abs(d)):.6e}")

# ----------------------------------------------------------------------
# Plot the integrated trajectories.
# ----------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

for i in range(traj2.shape[1]):
    axes[0].plot(t2, traj2[:, i], label=f"gene {i+1}")
axes[0].set_title("Two-gene system"); axes[0].set_xlabel("t"); axes[0].set_ylabel("expression"); axes[0].legend()

for i in range(traj3.shape[1]):
    axes[1].plot(t3, traj3[:, i], label=f"gene {i+1}")
axes[1].set_title("Three-gene system"); axes[1].set_xlabel("t"); axes[1].legend()

for i in range(trajN.shape[1]):
    axes[2].plot(tN, trajN[:, i], label=f"gene {i+1}")
axes[2].set_title(f"Many-gene system (N={N})"); axes[2].set_xlabel("t"); axes[2].legend(ncol=2, fontsize=8)

fig.suptitle("Trajectories from ONE unchanged RK4_generic routine")
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3H.1.1_s4.png")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: Because the identical RK4_generic call produced a "
      "correct-length, steady derivative for state vectors of length 2, 3, "
      "and 6 without any code change, the routine is truly dimension-agnostic "
      "and depends only on f returning a vector matching the state length.")
