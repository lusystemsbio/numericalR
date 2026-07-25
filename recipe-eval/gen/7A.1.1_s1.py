import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_diffusion(ngrid, X_all, dX, dt, D, P_init, nsteps):
    """
    Generic explicit forward-time centered-space (FTCS) solver for the
    1D diffusion equation dP/dt = D * d2P/dX2, no reaction term.

    Dirichlet boundaries: P = 0 at both ends (indices 0 and ngrid-1).
    Returns the distribution P advanced `nsteps` finite-difference steps.
    """
    # Stability factor (a.k.a. diffusion number / mesh Fourier number).
    # For the explicit FTCS scheme this must satisfy alpha <= 0.5.
    alpha = D * dt / dX**2

    # Work on a copy so the caller's initial condition is preserved.
    P = np.array(P_init, dtype=float)

    for _ in range(nsteps):
        P_next = P.copy()  # start from current state; boundaries stay fixed
        # Update every interior node from its two neighbors (centered space).
        for i in range(1, ngrid - 1):
            # Explicit forward-time update:
            # P_i_next = P_i + alpha*(P_{i+1} + P_{i-1} - 2*P_i)
            P_next[i] = P[i] + alpha * (P[i + 1] + P[i - 1] - 2.0 * P[i])
        # Enforce Dirichlet (P = 0) boundaries explicitly.
        P_next[0] = 0.0
        P_next[-1] = 0.0
        P = P_next

    return P, alpha


# ---------------------------------------------------------------------------
# Set up the spatial grid and parameters.
# ---------------------------------------------------------------------------
ngrid = 101
X_all = np.linspace(0.0, 1.0, ngrid)
dX = X_all[1] - X_all[0]
D = 1.0
# Choose dt so the stability factor stays small (well below 0.5).
dt = 0.2 * dX**2 / D
nsteps = 200

# Initial condition: a narrow peak in the middle, zero at the boundaries.
P0 = np.zeros(ngrid)
P0[ngrid // 2] = 1.0 / dX  # normalized-ish spike

# Advance the distribution several finite-difference steps.
P_final, alpha = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps)

print("Grid spacing dX:", dX)
print("Time step dt:", dt)
print("Diffusion coefficient D:", D)
print("Stability factor D*dt/dX^2 (must be <= 0.5):", alpha)
print("Number of finite-difference steps taken:", nsteps)
print("Initial total mass (sum(P)*dX):", np.sum(P0) * dX)
print("Final total mass (sum(P)*dX):", np.sum(P_final) * dX)
print("Final peak value P at center:", P_final[ngrid // 2])
print("Final min value of P:", np.min(P_final))
print("Final max value of P:", np.max(P_final))

# ---------------------------------------------------------------------------
# Separate check: the FTCS scheme is the diffusion analog of Euler integration.
# Advance P by ONE time step from its neighbors, two independent ways:
#   (1) the vectorized explicit update (Euler step on the discretized RHS),
#   (2) an explicit Euler step P_new = P + dt * (D * Laplacian(P)),
# where Laplacian(P)_i = (P_{i+1}+P_{i-1}-2P_i)/dX^2.
# If they agree, the finite-difference scheme is exactly Euler applied to
# the spatially discretized diffusion operator.
# ---------------------------------------------------------------------------
P_test = P0.copy()

# Method (1): one FTCS step using the stability factor form.
P_ftcs = P_test.copy()
P_ftcs[1:-1] = P_test[1:-1] + alpha * (P_test[2:] + P_test[:-2] - 2.0 * P_test[1:-1])
P_ftcs[0] = 0.0
P_ftcs[-1] = 0.0

# Method (2): one explicit Euler step on the discrete diffusion RHS.
laplacian = np.zeros(ngrid)
laplacian[1:-1] = (P_test[2:] + P_test[:-2] - 2.0 * P_test[1:-1]) / dX**2
P_euler = P_test + dt * (D * laplacian)
P_euler[0] = 0.0
P_euler[-1] = 0.0

max_diff = np.max(np.abs(P_ftcs - P_euler))
print("Max difference between one FTCS step and one Euler step:", max_diff)
print("FTCS == Euler (within 1e-12)?", bool(max_diff < 1e-12))
# Explanation: the two match because the FTCS update is literally
# P_new = P + dt * D * (discrete second derivative), i.e. a forward-Euler
# time step applied to the spatially discretized diffusion operator, which
# confirms the scheme is the diffusion analog of Euler ODE integration.

# ---------------------------------------------------------------------------
# Plot initial vs. advanced distribution.
# ---------------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(X_all, P0, label="P initial (t=0)", linestyle="--")
plt.plot(X_all, P_final, label=f"P after {nsteps} steps")
plt.xlabel("X")
plt.ylabel("P(X, t)")
plt.title("Explicit FTCS finite-difference solution of 1D diffusion")
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.1.1_s1.png")
