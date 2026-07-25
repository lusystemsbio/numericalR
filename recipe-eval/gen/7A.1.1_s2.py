import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_diffusion(ngrid, X_all, dX, dt, D, P_init, nsteps):
    """Explicit forward-time centered-space (FTCS) solver for dP/dt = D d2P/dX2.

    Dirichlet boundaries: P = 0 at both ends. Advances P by `nsteps`
    finite-difference steps and returns the final distribution plus the
    stability factor D*dt/dX^2 (must stay small for stability).
    """
    alpha = D * dt / dX**2  # stability factor
    P = P_init.copy().astype(float)
    for _ in range(nsteps):
        P_next = P.copy()  # start from current state (keeps Dirichlet ends = 0)
        # centered second difference over interior points only
        for i in range(1, ngrid - 1):
            # Euler-in-time update from the three neighbors P_{i-1}, P_i, P_{i+1}
            P_next[i] = P[i] + alpha * (P[i + 1] + P[i - 1] - 2.0 * P[i])
        P_next[0] = 0.0        # left Dirichlet boundary
        P_next[-1] = 0.0       # right Dirichlet boundary
        P = P_next
    return P, alpha


# ---- Grid and parameters --------------------------------------------------
ngrid = 101
X_all = np.linspace(0.0, 1.0, ngrid)
dX = X_all[1] - X_all[0]
D = 1.0
# choose dt so the stability factor D*dt/dX^2 stays small (< 0.5 required)
dt = 0.2 * dX**2 / D
nsteps = 200

# initial distribution: a narrow bump in the middle, zero at the boundaries
P0 = np.exp(-((X_all - 0.5) ** 2) / (2 * 0.05**2))
P0[0] = 0.0
P0[-1] = 0.0

# ---- Run the solver (exercised in 7A.2) -----------------------------------
P_final, alpha = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps)

print(f"Grid points ngrid: {ngrid}")
print(f"Grid spacing dX: {dX:.6f}")
print(f"Time step dt: {dt:.6e}")
print(f"Diffusion coefficient D: {D:.6f}")
print(f"Stability factor D*dt/dX^2 (must be < 0.5): {alpha:.6f}")
print(f"Number of finite-difference steps: {nsteps}")
print(f"Initial peak value P(mid, t=0): {P0[ngrid // 2]:.6f}")
print(f"Final peak value P(mid, t=nsteps*dt): {P_final[ngrid // 2]:.6f}")
print(f"Initial total mass (sum*dX): {np.sum(P0) * dX:.6f}")
print(f"Final total mass (sum*dX): {np.sum(P_final) * dX:.6f}")
print(f"Left boundary P[0] (should be 0): {P_final[0]:.6f}")
print(f"Right boundary P[-1] (should be 0): {P_final[-1]:.6f}")

# ---- Separate check: one explicit step == diffusion analog of Euler -------
# Advance a single time step "by hand" from the neighbors and compare against
# one call of the solver with nsteps=1. Agreement confirms the update is
# exactly Euler integration applied to the diffusion (centered-space) term.
i_check = ngrid // 2
manual_one_step = P0[i_check] + alpha * (
    P0[i_check + 1] + P0[i_check - 1] - 2.0 * P0[i_check]
)
P_one, _ = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, 1)
solver_one_step = P_one[i_check]

print(f"Manual single Euler step at mid: {manual_one_step:.10f}")
print(f"Solver single step at mid:       {solver_one_step:.10f}")
print(f"Difference (should be ~0): {abs(manual_one_step - solver_one_step):.3e}")
# One sentence: this matches because the FTCS update P_i + (D dt/dX^2)(P_{i+1}+P_{i-1}-2P_i)
# is literally forward-Euler time stepping (new = old + dt * rate) using the
# centered-space estimate of D*d2P/dX2 as the rate, so reproducing one step by
# hand confirms the scheme is the diffusion analog of Euler ODE integration.
print(
    "Check explanation: the scheme is new = old + dt*rate with rate = "
    "D*(P_{i+1}+P_{i-1}-2P_i)/dX^2, i.e. forward-Euler applied to diffusion, "
    "so the hand step matching the solver step confirms the result."
)

# ---- Plot -----------------------------------------------------------------
plt.figure(figsize=(8, 5))
plt.plot(X_all, P0, label="P at t = 0")
plt.plot(X_all, P_final, label=f"P after {nsteps} steps")
plt.xlabel("X")
plt.ylabel("P(X, t)")
plt.title("1D diffusion: explicit FTCS finite-difference scheme")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.1.1_s2.png")
