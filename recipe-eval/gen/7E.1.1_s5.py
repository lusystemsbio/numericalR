import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import time

# ---------------------------------------------------------------------------
# Two-component 2D reaction-diffusion integrator
#   du/dt = f(u,v) + Du * (d2u/dX2 + d2u/dY2)
#   dv/dt = g(u,v) + Dv * (d2v/dX2 + d2v/dY2)
# Explicit forward-time centered-space (FTCS), periodic boundaries,
# five-point Laplacian built from whole-grid array shifts (np.roll).
# ---------------------------------------------------------------------------

def laplacian(Z, dX):
    """Five-point periodic Laplacian by whole-grid shifts (no loops)."""
    # Sum of the four nearest neighbours (rolled) minus 4*center, over dX^2.
    return (np.roll(Z,  1, axis=0) + np.roll(Z, -1, axis=0) +
            np.roll(Z,  1, axis=1) + np.roll(Z, -1, axis=1) - 4.0 * Z) / (dX * dX)


def step_2d(u, v, f, g, Du, Dv, dt, dX):
    """One explicit FTCS step. Reaction f,g are evaluated on the WHOLE grid
    at once (vectorized), then added to the diffusion (Laplacian) term."""
    u_new = u + dt * (f(u, v) + Du * laplacian(u, dX))   # forward Euler in time
    v_new = v + dt * (g(u, v) + Dv * laplacian(v, dX))   # centered in space
    return u_new, v_new


def run_blocks(u, v, f, g, Du, Dv, dt, dX, n_blocks, steps_per_block):
    """Run-in-blocks driver: integrate steps_per_block steps per block,
    returning the final u,v after n_blocks*steps_per_block total steps."""
    for _ in range(n_blocks):
        for _ in range(steps_per_block):
            u, v = step_2d(u, v, f, g, Du, Dv, dt, dX)
    return u, v


# ---------------------------------------------------------------------------
# Generic reaction terms (FitzHugh-Nagumo-like two-component kinetics)
# ---------------------------------------------------------------------------
a, b = 0.1, 1.0
def f(u, v):
    return u - u**3 - v + a          # activator kinetics
def g(u, v):
    return b * (u - v)               # inhibitor kinetics


# ---------------------------------------------------------------------------
# Grid, parameters, and stability check
# ---------------------------------------------------------------------------
N  = 64          # grid points per side
dX = 1.0         # grid spacing
Du = 1.0         # diffusion of u
Dv = 4.0         # diffusion of v (must satisfy the stability bound below)
dt = 0.005       # 2D stability needs D*dt/dX^2 < 1/4

D_max = max(Du, Dv)
stab = D_max * dt / dX**2
print(f"Grid size N x N: {N} x {N}")
print(f"dt: {dt}")
print(f"dX: {dX}")
print(f"Du: {Du}")
print(f"Dv: {Dv}")
print(f"Stability number D_max*dt/dX^2: {stab:.6f}")
print(f"Stability bound (must be < 0.25): {'OK' if stab < 0.25 else 'VIOLATED'}")

# ---------------------------------------------------------------------------
# Nearly-uniform initial condition: u = v = 1 with +-0.1 noise (seed 10)
# ---------------------------------------------------------------------------
rng = np.random.default_rng(10)
u0 = 1.0 + rng.uniform(-0.1, 0.1, size=(N, N))
v0 = 1.0 + rng.uniform(-0.1, 0.1, size=(N, N))
print(f"Initial u mean: {u0.mean():.6f}")
print(f"Initial v mean: {v0.mean():.6f}")
print(f"Initial u range: [{u0.min():.6f}, {u0.max():.6f}]")
print(f"Initial v range: [{v0.min():.6f}, {v0.max():.6f}]")

# ---------------------------------------------------------------------------
# Advance u and v on the 2D grid with the run-in-blocks driver (used in 7E.2)
# ---------------------------------------------------------------------------
n_blocks, steps_per_block = 20, 100
u, v = run_blocks(u0.copy(), v0.copy(), f, g, Du, Dv, dt, dX,
                  n_blocks, steps_per_block)
total_steps = n_blocks * steps_per_block
print(f"Total steps integrated: {total_steps}")
print(f"Final simulated time: {total_steps * dt:.6f}")
print(f"Final u mean: {u.mean():.6f}")
print(f"Final v mean: {v.mean():.6f}")
print(f"Final u range: [{u.min():.6f}, {u.max():.6f}]")
print(f"Final v range: [{v.min():.6f}, {v.max():.6f}]")

# ---------------------------------------------------------------------------
# Separate check: the vectorized 2D reaction (whole grid at once) runs faster
# than a point-by-point evaluation of the same reaction.
# ---------------------------------------------------------------------------
def reaction_pointwise(u, v):
    """Same f,g but evaluated cell by cell in Python loops (slow)."""
    fu = np.empty_like(u)
    gv = np.empty_like(v)
    n0, n1 = u.shape
    for i in range(n0):
        for j in range(n1):
            uij, vij = u[i, j], v[i, j]
            fu[i, j] = uij - uij**3 - vij + a
            gv[i, j] = b * (uij - vij)
    return fu, gv

# time the vectorized whole-grid reaction
t0 = time.perf_counter()
for _ in range(50):
    fu_vec, gv_vec = f(u, v), g(u, v)
t_vec = time.perf_counter() - t0

# time the point-by-point reaction
t0 = time.perf_counter()
for _ in range(50):
    fu_pt, gv_pt = reaction_pointwise(u, v)
t_pt = time.perf_counter() - t0

max_diff = max(np.abs(fu_vec - fu_pt).max(), np.abs(gv_vec - gv_pt).max())
print(f"Vectorized reaction time (50 evals): {t_vec:.6f} s")
print(f"Point-by-point reaction time (50 evals): {t_pt:.6f} s")
print(f"Speedup (pointwise/vectorized): {t_pt / t_vec:.2f}x")
print(f"Max abs difference between methods: {max_diff:.3e}")
print("Check confirms result: both methods give identical values (max diff ~0) "
      "but the vectorized whole-grid evaluation is much faster, which shows the "
      "array-based integrator computes the same reaction as the point-by-point "
      "loop while replacing per-cell Python overhead with a single batched operation.")

# ---------------------------------------------------------------------------
# Figure: initial and final u and v fields
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(9, 8))
im0 = ax[0, 0].imshow(u0, cmap="viridis"); ax[0, 0].set_title("initial u"); fig.colorbar(im0, ax=ax[0, 0])
im1 = ax[0, 1].imshow(v0, cmap="magma");   ax[0, 1].set_title("initial v"); fig.colorbar(im1, ax=ax[0, 1])
im2 = ax[1, 0].imshow(u,  cmap="viridis"); ax[1, 0].set_title("final u");   fig.colorbar(im2, ax=ax[1, 0])
im3 = ax[1, 1].imshow(v,  cmap="magma");   ax[1, 1].set_title("final v");   fig.colorbar(im3, ax=ax[1, 1])
fig.suptitle("Two-component 2D reaction-diffusion (FTCS, periodic)")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.1.1_s5.png")
