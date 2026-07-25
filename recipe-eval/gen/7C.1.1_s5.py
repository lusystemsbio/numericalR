import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D, u0,
                                     nsteps=20000, record_every=None):
    """
    Generic n-component 1D reaction-diffusion integrator.

        du_k/dt = f_k(u) + D_k * d2u_k/dX2 ,   k = 0..n-1

    Vectorized explicit forward-time centered-space (FTCS) scheme with
    periodic (wrap-around) boundaries. All components advance together.

    Parameters
    ----------
    derivs : callable u -> (n, ngrid) array of reaction terms f_k(u)
    n      : number of coupled components
    ngrid  : number of spatial grid points
    dX     : spatial step
    dt     : time step
    D      : length-n vector of diffusion constants
    u0     : (n, ngrid) initial state
    """
    D = np.asarray(D, dtype=float).reshape(n, 1)   # broadcast over grid
    u = np.array(u0, dtype=float)                  # working copy, shape (n, ngrid)
    assert u.shape == (n, ngrid)

    inv_dX2 = 1.0 / (dX * dX)

    for step in range(nsteps):
        # Reaction term for every component at once -> (n, ngrid)
        react = derivs(u)

        # Centered-space Laplacian with wrap-around neighbors.
        # np.roll shifts along the grid axis, giving periodic boundaries.
        lap = (np.roll(u, -1, axis=1) - 2.0 * u + np.roll(u, 1, axis=1)) * inv_dX2

        # Explicit Euler update, each component with its own D_k.
        u = u + dt * (react + D * lap)

    return u


def dominant_wavelength(field, dX):
    """Wavelength (in X units) of the strongest nonzero spatial mode."""
    f = field - field.mean()
    spec = np.abs(np.fft.rfft(f))
    spec[0] = 0.0                                  # drop the DC component
    k = np.argmax(spec)                            # index of dominant mode
    ngrid = field.size
    return ngrid * dX / k                          # domain length / mode number


# --- Test model: two-component Turing (Schnakenberg) activator-inhibitor ---
# f_u = a - u + u^2 v ,  f_v = b - u^2 v
a, b = 0.1, 0.9


def schnakenberg(u):
    U, V = u[0], u[1]
    f_U = a - U + U * U * V
    f_V = b - U * U * V
    return np.stack([f_U, f_V])


n = 2
ngrid = 200
dX = 1.0
dt = 0.01
D = [1.0, 40.0]          # inhibitor (V) diffuses much faster -> Turing patterns
nsteps = 40000

# Homogeneous steady state: U* = a+b, V* = b/(a+b)^2
U_star = a + b
V_star = b / (a + b) ** 2
print(f"Homogeneous steady state U*: {U_star}")
print(f"Homogeneous steady state V*: {V_star}")

# Random initial perturbations about the steady state (sets pattern PHASE only).
rng = np.random.default_rng(12345)
u0 = np.empty((n, ngrid))
u0[0] = U_star + 0.01 * rng.standard_normal(ngrid)
u0[1] = V_star + 0.01 * rng.standard_normal(ngrid)

# Advance all components together on the wrap-around grid.
u_final = pde_fd_reaction_diffusion_multi(schnakenberg, n, ngrid, dX, dt, D,
                                          u0, nsteps=nsteps)

print(f"Number of components advanced: {n}")
print(f"Grid points: {ngrid}")
print(f"Steps integrated: {nsteps}")
print(f"Final U min: {u_final[0].min()}")
print(f"Final U max: {u_final[0].max()}")
print(f"Final V min: {u_final[1].min()}")
print(f"Final V max: {u_final[1].max()}")

lam_U = dominant_wavelength(u_final[0], dX)
lam_V = dominant_wavelength(u_final[1], dX)
print(f"Dominant wavelength of U pattern: {lam_U}")
print(f"Dominant wavelength of V pattern: {lam_V}")

# Wrap-around continuity check: neighbor difference across the periodic seam
seam_gap = abs(u_final[0, 0] - u_final[0, -1])
typical_gap = np.abs(np.diff(u_final[0])).mean()
print(f"Seam (wrap-around) gap in U: {seam_gap}")
print(f"Typical neighbor gap in U: {typical_gap}")

X = np.arange(ngrid) * dX
plt.figure(figsize=(8, 4))
plt.plot(X, u_final[0], label="U (activator)")
plt.plot(X, u_final[1], label="V (inhibitor)")
plt.xlabel("X")
plt.ylabel("concentration")
plt.title("n-component reaction-diffusion pattern (periodic grid)")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.1.1_s5.png")

# Why the check confirms the result:
print("Explanation: The pattern wavelength is fixed by the deterministic "
      "reaction-diffusion dispersion relation (the fastest-growing Turing mode), "
      "while only the pattern's position/phase depends on the random seed, so R "
      "and Python matching in wavelength but not position confirms both integrate "
      "the same coupled dynamics correctly.")
