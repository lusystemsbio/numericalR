import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Generic n-component 1D reaction-diffusion integrator (explicit FTCS)
#   du_k/dt = f_k(u) + D_k * d2u_k/dX2
# State is an (n, ngrid) matrix; D is a length-n vector; periodic boundaries.
# -----------------------------------------------------------------------------
def pde_fd_reaction_diffusion_multi(derivs, n, ngrid, dX, dt, D,
                                    u0, nsteps, record_every=None):
    """
    derivs : callable(u) -> (n, ngrid) array of reaction terms f_k(u)
    n      : number of coupled components
    ngrid  : number of spatial grid points
    dX     : spatial step
    dt     : time step
    D      : length-n vector of diffusion constants
    u0     : (n, ngrid) initial state matrix
    nsteps : number of forward-time steps to take
    Returns (u_final, times, snapshots) where snapshots is a list of recorded states.
    """
    D = np.asarray(D, dtype=float).reshape(n, 1)   # column vector broadcasts over grid
    u = np.array(u0, dtype=float).copy()           # working copy of the (n, ngrid) state
    inv_dX2 = 1.0 / (dX * dX)

    snapshots = []
    times = []
    if record_every is None:
        record_every = max(1, nsteps // 200)

    for step in range(nsteps):
        # --- Centered-space Laplacian with periodic (wrap-around) boundaries ---
        # np.roll wraps the last column to the front and vice-versa, giving
        # u[:, i+1] and u[:, i-1] with periodic edges, all components at once.
        u_right = np.roll(u, -1, axis=1)           # neighbour i+1 (wraps at right edge)
        u_left  = np.roll(u,  1, axis=1)           # neighbour i-1 (wraps at left edge)
        lap = (u_right - 2.0 * u + u_left) * inv_dX2   # (n, ngrid) second derivative

        # --- Reaction term for every component (vectorized user function) ---
        react = derivs(u)                          # (n, ngrid)

        # --- Forward-time explicit update, all components advanced together ---
        u = u + dt * (react + D * lap)

        if step % record_every == 0:
            snapshots.append(u.copy())
            times.append((step + 1) * dt)

    snapshots.append(u.copy())
    times.append(nsteps * dt)
    return u, np.array(times), snapshots


# -----------------------------------------------------------------------------
# Test system: Schnakenberg reaction kinetics (a classic Turing model).
#   f_u = a - u + u^2 v      (activator, small diffusion)
#   f_v = b - u^2 v          (inhibitor, large diffusion)
# Diffusion-driven instability => stationary pattern with a characteristic
# wavelength, seeded from random perturbations about the homogeneous state.
# -----------------------------------------------------------------------------
a, b = 0.1, 0.9

def schnakenberg(u):
    U = u[0]
    V = u[1]
    fU = a - U + U * U * V
    fV = b - U * U * V
    return np.vstack([fU, fV])

# Problem setup
n      = 2
ngrid  = 256
L      = 60.0
dX     = L / ngrid
dt     = 0.01
D      = np.array([1.0, 40.0])     # D_u << D_v drives the Turing instability
nsteps = 200000

# Homogeneous steady state: u* = a+b, v* = b/(a+b)^2
u_star = a + b
v_star = b / (u_star ** 2)
print(f"Homogeneous steady state u* = {u_star:.6f}")
print(f"Homogeneous steady state v* = {v_star:.6f}")

# Random initial perturbations about the steady state (fixed seed for this run)
rng = np.random.default_rng(12345)
u0 = np.empty((n, ngrid))
u0[0] = u_star + 0.01 * rng.standard_normal(ngrid)
u0[1] = v_star + 0.01 * rng.standard_normal(ngrid)

# CFL-type stability check for the explicit scheme (max over components)
cfl = D.max() * dt / (dX * dX)
print(f"Explicit-scheme stability number D_max*dt/dX^2 = {cfl:.6f} (must be < 0.5)")

# Integrate all components together on the wrap-around grid
u_final, times, snapshots = pde_fd_reaction_diffusion_multi(
    schnakenberg, n, ngrid, dX, dt, D, u0, nsteps)

# --- Report final-state summary numbers ---
print(f"Final activator U min = {u_final[0].min():.6f}")
print(f"Final activator U max = {u_final[0].max():.6f}")
print(f"Final activator U mean = {u_final[0].mean():.6f}")
print(f"Final inhibitor V min = {u_final[1].min():.6f}")
print(f"Final inhibitor V max = {u_final[1].max():.6f}")
print(f"Final inhibitor V mean = {u_final[1].mean():.6f}")

# --- Dominant wavelength of the emergent pattern (via FFT of activator) ---
# The pattern's wavelength is position-independent, so it is the quantity that
# should match between R and Python even when the random seeds place the
# stripes at different locations.
U = u_final[0] - u_final[0].mean()
spectrum = np.abs(np.fft.rfft(U))
freqs = np.fft.rfftfreq(ngrid, d=dX)     # cycles per unit length
kdom = np.argmax(spectrum[1:]) + 1       # skip the zero-frequency (mean) bin
dominant_wavelength = 1.0 / freqs[kdom]
n_stripes = L / dominant_wavelength
print(f"Dominant spatial mode index = {kdom}")
print(f"Dominant wavelength = {dominant_wavelength:.6f}")
print(f"Number of stripes across domain length L={L:.1f} = {n_stripes:.6f}")

# --- Plot ---
X = np.arange(ngrid) * dX
fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
ax[0].plot(X, u0[0], 'k--', alpha=0.5, label="U initial (random)")
ax[0].plot(X, u_final[0], 'b-', label="U final")
ax[0].set_ylabel("Activator U")
ax[0].legend(loc="upper right")
ax[0].set_title("Schnakenberg Turing pattern (n-component RD solver, periodic BC)")
ax[1].plot(X, u_final[1], 'r-', label="V final")
ax[1].set_ylabel("Inhibitor V")
ax[1].set_xlabel("X")
ax[1].legend(loc="upper right")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.1.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Check rationale: because the Turing wavelength is set by the reaction "
      "kinetics and diffusion constants (not by initial conditions), R and "
      "Python producing the same wavelength at different stripe positions "
      "confirms both integrators advance all coupled components identically on "
      "the wrap-around grid, with only the random seed shifting the phase.")
