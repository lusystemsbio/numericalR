import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps):
    """Explicit forward-time centered-space (FTCS) solver for dP/dt = D d2P/dX2.

    Dirichlet boundaries: P = 0 at both ends. Advances P by nsteps FD steps.
    """
    # Stability factor (must stay small, <= 0.5 for 1D explicit diffusion).
    r = D * dt / dX**2
    P = P0.astype(float).copy()
    for _ in range(nsteps):
        P_next = P.copy()  # start from current values
        # Update only interior points; boundaries stay clamped at 0 (Dirichlet).
        for i in range(1, ngrid - 1):
            # Explicit FTCS update: new value from current + curvature term.
            P_next[i] = P[i] + r * (P[i + 1] + P[i - 1] - 2 * P[i])
        # Enforce Dirichlet boundary conditions.
        P_next[0] = 0.0
        P_next[-1] = 0.0
        P = P_next
    return P, r


# --- Grid and parameters ---------------------------------------------------
ngrid = 101
X_all = np.linspace(0.0, 1.0, ngrid)
dX = X_all[1] - X_all[0]
D = 1.0
dt = 0.4 * dX**2 / D   # choose dt so stability factor r = 0.4 < 0.5
nsteps = 200

# Initial distribution: a localized bump in the middle, zero at boundaries.
P0 = np.exp(-((X_all - 0.5) ** 2) / (2 * 0.05**2))
P0[0] = 0.0
P0[-1] = 0.0

stability_factor = D * dt / dX**2
print(f"Number of grid points (ngrid): {ngrid}")
print(f"Grid spacing dX: {dX}")
print(f"Time step dt: {dt}")
print(f"Diffusion coefficient D: {D}")
print(f"Stability factor D*dt/dX^2: {stability_factor}")
print(f"Number of FD steps: {nsteps}")

# --- Advance the distribution several FD steps -----------------------------
P_final, r = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps)

print(f"Initial total mass (sum*dX): {np.sum(P0) * dX}")
print(f"Final total mass (sum*dX): {np.sum(P_final) * dX}")
print(f"Initial peak value: {np.max(P0)}")
print(f"Final peak value: {np.max(P_final)}")
print(f"Final peak location X: {X_all[np.argmax(P_final)]}")

# --- Separate check: one step is the diffusion analog of Euler integration --
# Euler for an ODE: y_next = y + dt * f(y). Here f(P)_i = D*(P_{i+1}+P_{i-1}-2P_i)/dX^2,
# the discrete Laplacian, so one FD step should equal one explicit Euler step.
i = 50  # an interior test point
rhs = D * (P0[i + 1] + P0[i - 1] - 2 * P0[i]) / dX**2  # f(P) at point i
euler_one_step = P0[i] + dt * rhs                       # explicit Euler update
P_one_step, _ = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, 1)  # one FD step
fd_one_step = P_one_step[i]

print(f"Euler-form one-step value at i={i}: {euler_one_step}")
print(f"FD scheme one-step value at i={i}: {fd_one_step}")
print(f"Difference (Euler vs FD one step): {abs(euler_one_step - fd_one_step)}")
# This check confirms the result because advancing P one step from its neighbors
# with P_next = P + dt*(D*Laplacian) is exactly explicit Euler applied to the
# spatially-discretized diffusion ODE, so identical values prove the scheme is
# the diffusion analog of Euler integration.

# --- Plot ------------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(X_all, P0, label="P initial")
plt.plot(X_all, P_final, label=f"P after {nsteps} FD steps")
plt.xlabel("X")
plt.ylabel("P(X, t)")
plt.title("Explicit FTCS finite-difference diffusion")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.1.1_s4.png")
