import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import time


# -------------------------------------------------------------------------
# Generic reaction terms f(u,v) and g(u,v).
# Any two-component kinetics can be plugged in here; we use a simple
# FitzHugh-Nagumo-like pair so the fields actually evolve.
# -------------------------------------------------------------------------
def f(u, v):
    return u - u**3 - v + 0.05

def g(u, v):
    return 0.08 * (u - 0.8 * v + 0.1)


# -------------------------------------------------------------------------
# Five-point Laplacian on a periodic grid, computed by whole-grid shifts.
# np.roll wraps the array, giving periodic boundary conditions for free.
# lap = (up + down + left + right - 4*center) / dX^2
# -------------------------------------------------------------------------
def laplacian(a, dX):
    return (np.roll(a, 1, axis=0) + np.roll(a, -1, axis=0) +
            np.roll(a, 1, axis=1) + np.roll(a, -1, axis=1) -
            4.0 * a) / (dX * dX)


# -------------------------------------------------------------------------
# One explicit forward-time centered-space (FTCS) step for both fields.
# The reaction is evaluated on the WHOLE grid at once (vectorized).
# -------------------------------------------------------------------------
def ftcs_step(u, v, Du, Dv, dt, dX):
    # reaction evaluated over the entire grid in a single call each
    du = f(u, v) + Du * laplacian(u, dX)
    dv = g(u, v) + Dv * laplacian(v, dX)
    # explicit Euler update in time
    u_new = u + dt * du
    v_new = v + dt * dv
    return u_new, v_new


# -------------------------------------------------------------------------
# Run-in-blocks driver: advance the fields in chunks of steps.
# -------------------------------------------------------------------------
def run_blocks(u, v, Du, Dv, dt, dX, n_blocks, steps_per_block):
    for _ in range(n_blocks):
        for _ in range(steps_per_block):
            u, v = ftcs_step(u, v, Du, Dv, dt, dX)
    return u, v


# -------------------------------------------------------------------------
# Set up the grid and parameters.
# -------------------------------------------------------------------------
N = 100          # grid points per side
dX = 1.0         # spatial step
Du = 1.0         # diffusion coefficient for u
Dv = 0.5         # diffusion coefficient for v

# 2D stability requires D*dt/dX^2 < 1/4  ->  dt < dX^2/(4*Dmax)
Dmax = max(Du, Dv)
dt = 0.005
stability_number = Dmax * dt / (dX * dX)
print(f"Stability number D*dt/dX^2 (must be < 0.25): {stability_number}")
print(f"Stability satisfied: {stability_number < 0.25}")

# nearly uniform initial condition: u = v = 1 with +/- 0.1 noise, seed 10
rng = np.random.default_rng(10)
u0 = 1.0 + (rng.random((N, N)) - 0.5) * 0.2   # 0.2 span -> +/- 0.1
v0 = 1.0 + (rng.random((N, N)) - 0.5) * 0.2
print(f"Initial u mean: {u0.mean()}")
print(f"Initial u min:  {u0.min()}")
print(f"Initial u max:  {u0.max()}")
print(f"Initial v mean: {v0.mean()}")
print(f"Initial v min:  {v0.min()}")
print(f"Initial v max:  {v0.max()}")

# -------------------------------------------------------------------------
# Advance the fields with the run-in-blocks driver.
# -------------------------------------------------------------------------
n_blocks = 20
steps_per_block = 100
total_steps = n_blocks * steps_per_block
u, v = run_blocks(u0.copy(), v0.copy(), Du, Dv, dt, dX, n_blocks, steps_per_block)

print(f"Total steps taken: {total_steps}")
print(f"Final time reached: {total_steps * dt}")
print(f"Final u mean: {u.mean()}")
print(f"Final u min:  {u.min()}")
print(f"Final u max:  {u.max()}")
print(f"Final v mean: {v.mean()}")
print(f"Final v min:  {v.min()}")
print(f"Final v max:  {v.max()}")


# -------------------------------------------------------------------------
# Speed check: vectorized whole-grid reaction vs. point-by-point loops.
# Both compute exactly one FTCS step of the reaction+diffusion update.
# -------------------------------------------------------------------------
ua = u0.copy()
va = v0.copy()

# --- vectorized: reaction evaluated on the whole grid at once ---
t0 = time.perf_counter()
for _ in range(10):
    ua, va = ftcs_step(ua, va, Du, Dv, dt, dX)
t_vec = time.perf_counter() - t0

# --- point-by-point: 1D-style nested Python loops over every cell ---
ub = u0.copy()
vb = v0.copy()
t0 = time.perf_counter()
for _ in range(10):
    unew = np.empty_like(ub)
    vnew = np.empty_like(vb)
    for i in range(N):
        ip = (i + 1) % N
        im = (i - 1) % N
        for j in range(N):
            jp = (j + 1) % N
            jm = (j - 1) % N
            lap_u = (ub[ip, j] + ub[im, j] + ub[i, jp] + ub[i, jm] - 4.0 * ub[i, j]) / (dX * dX)
            lap_v = (vb[ip, j] + vb[im, j] + vb[i, jp] + vb[i, jm] - 4.0 * vb[i, j]) / (dX * dX)
            unew[i, j] = ub[i, j] + dt * (f(ub[i, j], vb[i, j]) + Du * lap_u)
            vnew[i, j] = vb[i, j] + dt * (g(ub[i, j], vb[i, j]) + Dv * lap_v)
    ub, vb = unew, vnew
t_loop = time.perf_counter() - t0

max_diff = np.max(np.abs(ua - ub))
print(f"Vectorized 2D time (10 steps): {t_vec} s")
print(f"Point-by-point 1D-style time (10 steps): {t_loop} s")
print(f"Speedup factor (loop / vectorized): {t_loop / t_vec}")
print(f"Max difference between the two methods: {max_diff}")
print(f"Vectorized is faster: {t_vec < t_loop}")
# Explanation: the two methods produce numerically identical fields
# (max difference ~ machine precision) yet the vectorized version is far
# faster, confirming that evaluating the reaction on the whole grid at once
# gives the same result as the point-by-point approach but with less overhead.

# -------------------------------------------------------------------------
# Plot final u and v fields.
# -------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 5))
im0 = axes[0].imshow(u, origin="lower", cmap="viridis")
axes[0].set_title(f"u at t = {total_steps * dt}")
fig.colorbar(im0, ax=axes[0])
im1 = axes[1].imshow(v, origin="lower", cmap="magma")
axes[1].set_title(f"v at t = {total_steps * dt}")
fig.colorbar(im1, ax=axes[1])
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.1.1_s2.png")
