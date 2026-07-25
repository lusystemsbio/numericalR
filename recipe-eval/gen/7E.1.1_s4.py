import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import time

# ---------------------------------------------------------------------------
# Generic two-component 2D reaction-diffusion system:
#   du/dt = f(u,v) + Du*(d2u/dX2 + d2u/dY2)
#   dv/dt = g(u,v) + Dv*(d2v/dX2 + d2v/dY2)
# Explicit forward-time centered-space (FTCS) integrator, periodic boundaries,
# five-point Laplacian via whole-grid array shifts.
# ---------------------------------------------------------------------------

# Reaction terms (generic; a simple activator-inhibitor style choice)
def f(u, v):
    return u - u**3 - v          # activator reaction
def g(u, v):
    return 0.25 * (u - v)        # inhibitor reaction

# Five-point Laplacian using np.roll shifts (periodic wrap), whole grid at once.
def laplacian(A, dX):
    # shift up/down (axis 0) and left/right (axis 1), sum neighbors, subtract 4*center
    return (np.roll(A, 1, axis=0) + np.roll(A, -1, axis=0) +
            np.roll(A, 1, axis=1) + np.roll(A, -1, axis=1) - 4.0 * A) / dX**2

# One explicit FTCS step, evaluating the reaction on the whole grid at once.
def step(u, v, Du, Dv, dt, dX):
    u_new = u + dt * (f(u, v) + Du * laplacian(u, dX))   # update u everywhere
    v_new = v + dt * (g(u, v) + Dv * laplacian(v, dX))   # update v everywhere
    return u_new, v_new

# Run-in-blocks driver: advance 'nsteps' total, reported in blocks.
def run_in_blocks(u, v, Du, Dv, dt, dX, nsteps, nblocks):
    per_block = nsteps // nblocks
    for b in range(nblocks):
        for _ in range(per_block):
            u, v = step(u, v, Du, Dv, dt, dX)
    return u, v

# ---------------------------------------------------------------------------
# Set up grid, parameters, and stability check
# ---------------------------------------------------------------------------
N = 64
dX = 1.0
Du, Dv = 0.5, 0.25
dt = 0.005                      # 2D stability: D*dt/dX^2 < 1/4
Dmax = max(Du, Dv)
stability = Dmax * dt / dX**2
print(f"Stability ratio D*dt/dX^2 = {stability:.6f} (must be < 0.25)")
print(f"Stable: {stability < 0.25}")

# Nearly uniform initial condition u = v = 1 with +-0.1 noise, seed 10
rng = np.random.default_rng(10)
u0 = 1.0 + 0.1 * (2 * rng.random((N, N)) - 1)
v0 = 1.0 + 0.1 * (2 * rng.random((N, N)) - 1)

print(f"Initial u mean = {u0.mean():.6f}, min = {u0.min():.6f}, max = {u0.max():.6f}")
print(f"Initial v mean = {v0.mean():.6f}, min = {v0.min():.6f}, max = {v0.max():.6f}")

# Advance the system (exercised further in 7E.2)
nsteps, nblocks = 2000, 10
u, v = run_in_blocks(u0.copy(), v0.copy(), Du, Dv, dt, dX, nsteps, nblocks)

print(f"After {nsteps} steps: u mean = {u.mean():.6f}, min = {u.min():.6f}, max = {u.max():.6f}")
print(f"After {nsteps} steps: v mean = {v.mean():.6f}, min = {v.min():.6f}, max = {v.max():.6f}")

# ---------------------------------------------------------------------------
# Separate check: vectorized whole-grid reaction vs point-by-point loop
# ---------------------------------------------------------------------------
# Vectorized: evaluate reaction on the entire grid in one array operation.
t0 = time.perf_counter()
for _ in range(50):
    _ = f(u, v)
t_vec = time.perf_counter() - t0

# Point-by-point: emulate a 1D-style approach looping over every grid point.
def f_scalar(uu, vv):
    return uu - uu**3 - vv
t0 = time.perf_counter()
for _ in range(50):
    out = np.empty_like(u)
    for i in range(N):
        for j in range(N):
            out[i, j] = f_scalar(u[i, j], v[i, j])
t_loop = time.perf_counter() - t0

print(f"Vectorized reaction time (50 evals) = {t_vec:.6f} s")
print(f"Point-by-point reaction time (50 evals) = {t_loop:.6f} s")
print(f"Speedup factor (loop / vectorized) = {t_loop / t_vec:.2f}")
print(f"Vectorized faster: {t_vec < t_loop}")
# Explanation: this check confirms the result because both methods must produce
# the identical reaction field, so the vectorized version being much faster shows
# that evaluating the whole grid at once is a pure performance win with no change
# in the computed values.
print(f"Max abs difference between methods = {np.max(np.abs(f(u, v) - out)):.3e}")

# ---------------------------------------------------------------------------
# Plot final u and v fields
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
im0 = axes[0].imshow(u, cmap="viridis", origin="lower")
axes[0].set_title("u after integration")
fig.colorbar(im0, ax=axes[0])
im1 = axes[1].imshow(v, cmap="magma", origin="lower")
axes[1].set_title("v after integration")
fig.colorbar(im1, ax=axes[1])
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.1.1_s4.png")
