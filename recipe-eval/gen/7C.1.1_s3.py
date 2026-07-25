import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D, u0, nsteps):
    """
    Generic n-component 1D reaction-diffusion integrator.

        du_k/dt = f_k(u) + D_k * d2u_k/dX2   for k = 0..n-1

    Vectorized explicit forward-time centered-space (FTCS) scheme with
    periodic (wrap-around) boundaries. Advances ALL components together.

    Parameters
    ----------
    derivs : callable(u) -> (n, ngrid) array of reaction terms f_k(u)
    n      : number of components
    ngrid  : number of spatial grid points
    dX     : spatial step
    dt     : time step
    D      : length-n vector of diffusion constants
    u0     : (n, ngrid) initial state matrix
    nsteps : number of time steps to advance
    """
    D = np.asarray(D, dtype=float).reshape(n, 1)   # column so it broadcasts over grid
    u = np.array(u0, dtype=float)                  # working copy, shape (n, ngrid)

    for _ in range(nsteps):
        # Centered second difference with wrap-around neighbours.
        # np.roll shifts along the grid axis (axis=1) periodically.
        lap = (np.roll(u, -1, axis=1) - 2.0 * u + np.roll(u, +1, axis=1)) / dX**2

        # Reaction term for every component at once, shape (n, ngrid).
        react = derivs(u)

        # Explicit Euler update: all components advanced simultaneously.
        u = u + dt * (react + D * lap)

    return u


# --- Test: 2-component activator-inhibitor (Turing) system ---------------
# Linear reaction kinetics chosen so a band of wavenumbers is unstable,
# giving a pattern with a well-defined preferred wavelength.
a, b, c, d = 1.0, -1.5, 2.0, -1.0   # Jacobian entries of the reaction

def turing_derivs(u):
    U, V = u[0], u[1]
    fU = a * U + b * V
    fV = c * U + d * V
    return np.vstack([fU, fV])

# Grid and integration parameters.
n = 2
ngrid = 200
L = 40.0
dX = L / ngrid
D = np.array([1.0, 20.0])           # activator slow, inhibitor fast -> Turing patterns
# Stable explicit step: dt < dX^2 / (2*max(D)).
dt = 0.2 * dX**2 / (2.0 * D.max())
nsteps = 40000

# Random initial perturbations about the (0,0) homogeneous state.
rng = np.random.default_rng(12345)
u0 = 0.01 * rng.standard_normal((n, ngrid))

uT = pde_fd_reaction_diffusion_multi(turing_derivs, n, ngrid, dX, dt, D, u0, nsteps)

# --- Analyze the emergent pattern wavelength via the FFT ------------------
Uf = uT[0] - uT[0].mean()
spec = np.abs(np.fft.rfft(Uf))
freqs = np.fft.rfftfreq(ngrid, d=dX)     # cycles per unit length
kmax = np.argmax(spec[1:]) + 1           # skip DC component
dom_freq = freqs[kmax]
dom_wavelength = 1.0 / dom_freq
n_modes = dom_freq * L                   # number of full waves around the ring

print(f"grid points ngrid: {ngrid}")
print(f"domain length L: {L}")
print(f"spatial step dX: {dX}")
print(f"time step dt: {dt}")
print(f"number of steps nsteps: {nsteps}")
print(f"diffusion constants D: {D.tolist()}")
print(f"final min of component 0 (U): {uT[0].min()}")
print(f"final max of component 0 (U): {uT[0].max()}")
print(f"final min of component 1 (V): {uT[1].min()}")
print(f"final max of component 1 (V): {uT[1].max()}")
print(f"dominant spatial frequency (cycles/length): {dom_freq}")
print(f"dominant wavelength: {dom_wavelength}")
print(f"number of full waves around ring: {n_modes}")

# --- Plot the two components on the wrap-around grid ----------------------
X = np.arange(ngrid) * dX
plt.figure(figsize=(9, 4))
plt.plot(X, uT[0], label="component 0 (activator U)")
plt.plot(X, uT[1], label="component 1 (inhibitor V)")
plt.xlabel("X")
plt.ylabel("u_k")
plt.title(f"n-component reaction-diffusion (periodic), wavelength ~ {dom_wavelength:.2f}")
plt.legend()
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.1.1_s3.png")

# --- Why this check confirms the result ----------------------------------
print("Explanation: Because the pattern's wavelength is set by the "
      "deterministic Turing instability (the diffusion constants and reaction "
      "Jacobian) while its spatial phase is set by the random seed, R and "
      "Python producing the same wavelength at different positions confirms "
      "the solver advances all components together correctly and is not "
      "seed-dependent in its physics.")
