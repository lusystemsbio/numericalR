import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps):
    """Explicit forward-time centered-space (FTCS) solver for dP/dt = D d2P/dX2.

    Advances the initial distribution P0 by nsteps finite-difference steps,
    with Dirichlet boundaries P = 0 held fixed at both ends.
    """
    P = P0.copy().astype(float)
    P[0] = 0.0            # Dirichlet boundary at left end
    P[-1] = 0.0           # Dirichlet boundary at right end
    alpha = D * dt / dX**2   # stability factor; must stay small (<= 0.5)

    for _ in range(nsteps):
        P_next = P.copy()
        # update only interior points; boundaries stay pinned at 0
        for i in range(1, ngrid - 1):
            # FTCS explicit update: new value from current value + neighbors
            P_next[i] = P[i] + alpha * (P[i + 1] + P[i - 1] - 2.0 * P[i])
        P = P_next
    return P


# ---- Problem setup ----
ngrid = 101
X_all = np.linspace(0.0, 1.0, ngrid)
dX = X_all[1] - X_all[0]
D = 1.0
dt = 0.4 * dX**2 / D          # choose dt so stability factor = 0.4 < 0.5
alpha = D * dt / dX**2

print(f"Number of grid points ngrid: {ngrid}")
print(f"Grid spacing dX: {dX:.6f}")
print(f"Time step dt: {dt:.8f}")
print(f"Diffusion coefficient D: {D}")
print(f"Stability factor D*dt/dX^2 (alpha): {alpha:.6f}")
print(f"Stability satisfied (alpha <= 0.5): {alpha <= 0.5}")

# ---- Initial distribution: a localized pulse in the interior ----
P0 = np.zeros(ngrid)
center = ngrid // 2
P0[center - 5:center + 6] = 1.0   # box pulse
P0[0] = 0.0
P0[-1] = 0.0

print(f"Initial total mass (sum P0 * dX): {np.sum(P0) * dX:.6f}")
print(f"Initial peak value: {np.max(P0):.6f}")

# ---- Advance several finite-difference steps (as exercised in 7A.2) ----
nsteps = 200
P = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, nsteps)

print(f"Number of finite-difference steps advanced: {nsteps}")
print(f"Total mass after diffusion (sum P * dX): {np.sum(P) * dX:.6f}")
print(f"Peak value after diffusion: {np.max(P):.6f}")
print(f"Peak location X after diffusion: {X_all[np.argmax(P)]:.6f}")
print(f"P at left boundary (should be 0): {P[0]:.6e}")
print(f"P at right boundary (should be 0): {P[-1]:.6e}")

# ---- Separate check: the FTCS update is the diffusion analog of Euler ----
# Euler for an ODE: y_next = y + dt * f(y).  Here f = D * d2P/dX2 approximated
# by the centered second difference, so one manual step should equal one
# solver step from the same state.
i = center                          # test one interior point
d2P = (P0[i + 1] + P0[i - 1] - 2.0 * P0[i]) / dX**2   # centered 2nd difference
euler_step = P0[i] + dt * D * d2P                     # explicit Euler form
solver_one_step = pde_fd_diffusion(ngrid, X_all, dX, dt, D, P0, 1)[i]

print(f"Euler-form one-step value at center: {euler_step:.10f}")
print(f"Solver one-step value at center: {solver_one_step:.10f}")
print(f"Euler-vs-solver difference at center: {abs(euler_step - solver_one_step):.3e}")
# One sentence: this matches because writing the FTCS update as
# P + dt*(D*d2P/dX2) is exactly Euler's y_next = y + dt*f(y) with f being the
# diffusion operator, so agreement confirms the scheme advances P one Euler
# time step from its neighbors.
print("Check: FTCS equals Euler integration of the diffusion operator ->",
      np.isclose(euler_step, solver_one_step))

# ---- Plot ----
plt.figure(figsize=(8, 5))
plt.plot(X_all, P0, label="initial P(X, 0)", linestyle="--")
plt.plot(X_all, P, label=f"P(X) after {nsteps} steps")
plt.xlabel("X")
plt.ylabel("P")
plt.title("Explicit FTCS diffusion (Dirichlet boundaries)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7A.1.1_s5.png")
