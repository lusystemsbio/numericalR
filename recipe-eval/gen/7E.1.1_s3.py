import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import time


# --- Reaction terms for a generic two-component system ---------------------
# We use a Gierer-Meinhardt-like activator-inhibitor kinetics as a concrete
# instance of the generic f(u,v), g(u,v).  The integrator itself is generic:
# any f, g may be substituted.
def f(u, v):
    # activator reaction
    return u * u / v - u


def g(u, v):
    # inhibitor reaction
    return u * u - v


# --- Five-point periodic Laplacian via whole-grid array shifts -------------
def laplacian(a, dX):
    # np.roll implements periodic boundaries by wrapping the grid.
    lap = (np.roll(a, 1, axis=0) + np.roll(a, -1, axis=0) +
           np.roll(a, 1, axis=1) + np.roll(a, -1, axis=1) - 4.0 * a)
    return lap / (dX * dX)


# --- One explicit forward-time centered-space (FTCS) step ------------------
def step(u, v, Du, Dv, dt, dX):
    # reaction evaluated on the WHOLE grid at once (vectorized)
    u_new = u + dt * (f(u, v) + Du * laplacian(u, dX))
    v_new = v + dt * (g(u, v) + Dv * laplacian(v, dX))
    return u_new, v_new


# --- Run-in-blocks driver --------------------------------------------------
def integrate(u, v, Du, Dv, dt, dX, n_blocks, steps_per_block):
    for _ in range(n_blocks):
        for _ in range(steps_per_block):
            u, v = step(u, v, Du, Dv, dt, dX)
    return u, v


# --- Set up the generic solver test ----------------------------------------
N = 64                 # grid points per side
dX = 1.0               # grid spacing
Du = 1.0               # diffusion coefficient of u
Dv = 4.0               # diffusion coefficient of v

# 2D stability requires D*dt/dX^2 < 1/4.
dt = 0.005
Dmax = max(Du, Dv)
stability_number = Dmax * dt / (dX * dX)
print("Stability number D*dt/dX^2 =", stability_number)
print("Stability satisfied (< 0.25):", stability_number < 0.25)

# Nearly uniform initial condition u = v = 1 with +-0.1 noise, seed 10.
rng = np.random.default_rng(10)
u = 1.0 + 0.2 * (rng.random((N, N)) - 0.5)   # uniform in [-0.1, +0.1] about 1
v = 1.0 + 0.2 * (rng.random((N, N)) - 0.5)

print("Initial u mean:", u.mean())
print("Initial v mean:", v.mean())
print("Initial u min/max:", u.min(), u.max())
print("Initial v min/max:", v.min(), v.max())

# Advance the fields on the 2D grid using the run-in-blocks driver.
n_blocks = 20
steps_per_block = 50
u, v = integrate(u, v, Du, Dv, dt, dX, n_blocks, steps_per_block)
total_steps = n_blocks * steps_per_block

print("Total time steps advanced:", total_steps)
print("Final simulated time:", total_steps * dt)
print("Final u mean:", u.mean())
print("Final v mean:", v.mean())
print("Final u min/max:", u.min(), u.max())
print("Final v min/max:", v.min(), v.max())


# --- Performance check: vectorized whole-grid vs point-by-point ------------
# The vectorized integrator evaluates the reaction on the entire grid in one
# array operation; the point-by-point approach loops over every cell.  Both
# compute exactly the same reaction values, so timing them on one reaction
# evaluation is a fair comparison of the two evaluation strategies.
def reaction_vectorized(u, v):
    return f(u, v), g(u, v)


def reaction_pointwise(u, v):
    # 1D-style scalar loop over each grid cell
    fu = np.empty_like(u)
    gv = np.empty_like(v)
    ny, nx = u.shape
    for i in range(ny):
        for j in range(nx):
            ui = u[i, j]
            vi = v[i, j]
            fu[i, j] = ui * ui / vi - ui
            gv[i, j] = ui * ui - vi
    return fu, gv

t0 = time.perf_counter()
fu_vec, gv_vec = reaction_vectorized(u, v)
t1 = time.perf_counter()
fu_pt, gv_pt = reaction_pointwise(u, v)
t2 = time.perf_counter()

vec_time = t1 - t0
pt_time = t2 - t1
print("Vectorized reaction time (s):", vec_time)
print("Point-by-point reaction time (s):", pt_time)
print("Speedup (pointwise / vectorized):", pt_time / vec_time)
print("Max abs difference between methods:", np.max(np.abs(fu_vec - fu_pt)))
print("Vectorized faster:", vec_time < pt_time)

# One-sentence explanation of why this check confirms the result:
# Because both methods produce identical reaction values (max difference ~0)
# yet the whole-grid array evaluation runs in far less time, the check
# confirms that vectorizing the reaction over the 2D grid gives the same
# physics as a per-point loop while being faster.
print("Why: identical outputs at much lower runtime confirm the vectorized "
      "whole-grid reaction is equivalent yet faster than the per-point loop.")


# --- Visualize the advanced 2D fields --------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
im0 = axes[0].imshow(u, cmap="viridis", origin="lower")
axes[0].set_title("u field (final)")
fig.colorbar(im0, ax=axes[0])
im1 = axes[1].imshow(v, cmap="magma", origin="lower")
axes[1].set_title("v field (final)")
fig.colorbar(im1, ax=axes[1])
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.1.1_s3.png")
