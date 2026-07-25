import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Generic explicit forward-time centered-space (FTCS) diffusion solver.
# Solves dP/dt = D * d2P/dX2 with Dirichlet boundaries P=0 at both ends.
# ----------------------------------------------------------------------
def pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps):
    # stability (diffusion) factor; must stay small (<= 0.5 for 1D FTCS)
    alpha = D * dt / dX**2

    # work on a copy so the caller's initial condition is preserved
    P = P0.astype(float).copy()

    for _ in range(nsteps):
        P_next = P.copy()  # start from current state
        # update interior nodes; boundaries stay fixed (Dirichlet P=0)
        for i in range(1, ngrid - 1):
            # explicit FTCS update: Euler-in-time, centered-in-space Laplacian
            P_next[i] = P[i] + alpha * (P[i + 1] + P[i - 1] - 2.0 * P[i])
        # enforce Dirichlet boundary conditions
        P_next[0] = 0.0
        P_next[-1] = 0.0
        P = P_next  # advance one full time step

    return P, alpha


# ----------------------------------------------------------------------
# Problem setup / grid
# ----------------------------------------------------------------------
ngrid = 101
X_all = np.linspace(0.0, 1.0, ngrid)
dX = X_all[1] - X_all[0]
D = 1.0
dt = 0.4 * dX**2 / D          # chosen so D*dt/dX^2 = 0.4 (< 0.5, stable)
nsteps = 200

# initial distribution: a localized bump in the interior, zero at the ends
P0 = np.exp(-((X_all - 0.5) ** 2) / (2 * 0.05**2))
P0[0] = 0.0
P0[-1] = 0.0

# ----------------------------------------------------------------------
# Run the solver for several finite-difference steps (used by 7A.2)
# ----------------------------------------------------------------------
P_final, alpha = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps)

print(f"Number of grid points ngrid: {ngrid}")
print(f"Grid spacing dX: {dX:.6f}")
print(f"Time step dt: {dt:.8f}")
print(f"Diffusion coefficient D: {D:.6f}")
print(f"Stability factor D*dt/dX^2 (alpha): {alpha:.6f}")
print(f"Stable (alpha <= 0.5): {alpha <= 0.5}")
print(f"Number of time steps advanced: {nsteps}")
print(f"Initial peak value P0 max: {P0.max():.6f}")
print(f"Final peak value P max after {nsteps} steps: {P_final.max():.6f}")
print(f"Sum of P initial (mass proxy): {P0.sum():.6f}")
print(f"Sum of P final (mass proxy): {P_final.sum():.6f}")
print(f"P_final at left boundary: {P_final[0]:.6e}")
print(f"P_final at right boundary: {P_final[-1]:.6e}")

# ----------------------------------------------------------------------
# Separate check: the scheme is the diffusion analog of Euler integration.
# Advance P by ONE step from its neighbors, at a single test node, by hand,
# and compare with one step of the full solver.
# ----------------------------------------------------------------------
i_test = ngrid // 2
# one explicit Euler-like step at node i_test: P_next = P + D*(dt/dX^2)*Laplacian
manual_one_step = P0[i_test] + D * (dt / dX**2) * (
    P0[i_test + 1] + P0[i_test - 1] - 2.0 * P0[i_test]
)
P_one_step, _ = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, 1)
solver_one_step = P_one_step[i_test]

print(f"Manual single-step update at node {i_test}: {manual_one_step:.10f}")
print(f"Solver single-step update at node {i_test}: {solver_one_step:.10f}")
print(f"Absolute difference (manual vs solver): {abs(manual_one_step - solver_one_step):.3e}")

# Explanation (one sentence):
# The manual one-step value matches the solver because P_i_next = P_i +
# D*(dt/dX^2)*(P_{i+1}+P_{i-1}-2*P_i) is exactly Euler's method dP/dt ~=
# (P_next - P)/dt applied with the centered-difference Laplacian, so
# reproducing one step by hand confirms the scheme advances P from its
# neighbors just like Euler integration for ODEs.
print("Check: single explicit step equals Euler integration using the "
      "centered-difference Laplacian of P at its neighbors.")

# ----------------------------------------------------------------------
# Plot initial vs final distribution
# ----------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(X_all, P0, label="Initial P (t=0)")
plt.plot(X_all, P_final, label=f"P after {nsteps} steps")
plt.xlabel("X")
plt.ylabel("P(X, t)")
plt.title("1D Diffusion: Explicit FTCS Finite-Difference Scheme")
plt.legend()
plt.grid(True)
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.1.1_s3.png")
