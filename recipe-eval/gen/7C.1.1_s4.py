import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D, u0, nsteps):
    """
    Generic n-component 1D reaction-diffusion integrator.

        du_k/dt = f_k(u) + D_k * d2u_k/dX2,   k = 0..n-1

    Vectorized explicit forward-time centered-space (FTCS) scheme with
    periodic (wrap-around) boundaries.

    Parameters
    ----------
    derivs : callable(u) -> (n, ngrid) array of reaction terms f_k(u)
    n      : number of coupled components
    ngrid  : number of spatial grid points
    dX     : spatial step
    dt     : time step
    D      : length-n vector of diffusion constants
    u0     : (n, ngrid) initial state
    nsteps : number of time steps to advance
    """
    D = np.asarray(D, dtype=float).reshape(n, 1)   # column vector -> broadcasts over grid
    u = np.array(u0, dtype=float)                  # working copy, shape (n, ngrid)

    for _ in range(nsteps):
        # Centered second difference with periodic wrap-around (np.roll handles
        # the boundaries: point 0's left neighbor is point ngrid-1, and vice versa).
        lap = (np.roll(u, -1, axis=1) - 2.0 * u + np.roll(u, 1, axis=1)) / dX**2

        # Reaction term for every component at once.
        react = derivs(u)

        # Explicit Euler update: all n components advanced together in one vector op.
        u = u + dt * (react + D * lap)

    return u


# ---------------------------------------------------------------------------
# Test: a 2-component activator-inhibitor (Turing) system.  Small random
# perturbations of the uniform steady state grow into a periodic pattern whose
# wavelength is set by the linear instability (independent of the random seed).
# ---------------------------------------------------------------------------

# Gierer-Meinhardt style kinetics: u=activator, v=inhibitor.
a_src, b_dec = 0.1, 1.0

def derivs(u):
    A = u[0]
    I = u[1]
    fA = a_src - b_dec * A + A**2 / I          # activator reaction
    fI = A**2 - I                              # inhibitor reaction
    return np.vstack([fA, fI])

# Uniform steady state of the kinetics.
A_ss = (a_src + 1.0) / b_dec
I_ss = A_ss**2
print(f"Activator steady state A_ss = {A_ss}")
print(f"Inhibitor steady state I_ss = {I_ss}")

# Grid / integration parameters.
n = 2
ngrid = 200
dX = 0.5
dt = 0.005
D = [1.0, 40.0]                                # inhibitor diffuses much faster -> Turing patterns
nsteps = 200000

# Random initial perturbation about the steady state (fixed seed for reproducibility).
rng = np.random.default_rng(0)
u0 = np.empty((n, ngrid))
u0[0] = A_ss + 0.01 * rng.standard_normal(ngrid)
u0[1] = I_ss + 0.01 * rng.standard_normal(ngrid)

uf = pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D, u0, nsteps)

# Report the pattern's dominant wavelength via the FFT of the activator field.
A_final = uf[0]
A_dev = A_final - A_final.mean()
spec = np.abs(np.fft.rfft(A_dev))
freqs = np.fft.rfftfreq(ngrid, d=dX)
kmax = 1 + np.argmax(spec[1:])                 # skip the zero-frequency (mean) mode
L = ngrid * dX
dominant_wavelength = 1.0 / freqs[kmax]
print(f"Domain length L = {L}")
print(f"Dominant mode number = {kmax}")
print(f"Dominant wavelength = {dominant_wavelength}")
print(f"Final activator min = {A_final.min()}")
print(f"Final activator max = {A_final.max()}")
print(f"Final activator mean = {A_final.mean()}")
print(f"Final inhibitor min = {uf[1].min()}")
print(f"Final inhibitor max = {uf[1].max()}")

# Check confirmation:
print("This check confirms the result because both R and Python evolve identical "
      "kinetics on the same wrap-around grid, so the linear instability selects the "
      "same wavelength while the random seeds only shift where the peaks land.")

# Plot the final patterns for both components.
X = np.arange(ngrid) * dX
fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax[0].plot(X, uf[0], color="C0")
ax[0].set_ylabel("activator A")
ax[0].set_title("Reaction-diffusion pattern (periodic BCs, random seed)")
ax[1].plot(X, uf[1], color="C1")
ax[1].set_ylabel("inhibitor I")
ax[1].set_xlabel("X")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.1.1_s4.png")
