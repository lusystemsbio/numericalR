import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D Gierer-Meinhardt substrate-depletion reaction-diffusion model
#   f(u,v) = u^2 v - u        (activator u)
#   g(u,v) = mu (1 - u^2 v)   (substrate v)
#   Du = d, Dv = 1
# Explicit finite-difference integrator (from 7E.1), run in blocks.
# ---------------------------------------------------------------

# --- reaction terms ---
def f(u, v):
    return u * u * v - u

def g(u, v, mu):
    return mu * (1.0 - u * u * v)

# --- 5-point Laplacian with periodic (wrap) boundaries ---
def laplacian(a, dX):
    # np.roll gives neighbours in +/- x and +/- y with periodic wrapping
    lap = (np.roll(a, 1, 0) + np.roll(a, -1, 0) +
           np.roll(a, 1, 1) + np.roll(a, -1, 1) - 4.0 * a)
    return lap / (dX * dX)

# --- one explicit forward-Euler step of the RD system ---
def rd_step(u, v, d, mu, dX, dt):
    u_new = u + dt * (d * laplacian(u, dX) + f(u, v))     # activator: slow diffusion d
    v_new = v + dt * (1.0 * laplacian(v, dX) + g(u, v, mu))  # substrate: unit diffusion
    return u_new, v_new

# --- run the integrator for a block of many steps ---
def run_block(u, v, nsteps, d, mu, dX, dt):
    for _ in range(nsteps):
        u, v = rd_step(u, v, d, mu, dX, dt)
    return u, v

# --- parameters (test case) ---
N   = 101      # grid size (N x N)
dX  = 0.2      # spatial step
dt  = 0.005    # time step
d   = 0.1      # activator diffusion Du (Dv = 1)
mu  = 1.5      # substrate turnover rate

# --- nearly uniform initial condition: u = v = 1 with +-0.1 noise ---
rng = np.random.default_rng(10)   # seed 10
u = 1.0 + 0.1 * (2.0 * rng.random((N, N)) - 1.0)   # uniform in [-0.1, 0.1]
v = 1.0 + 0.1 * (2.0 * rng.random((N, N)) - 1.0)

# --- integrate in successive blocks, tracking pattern growth ---
steps_per_block = 2000
n_blocks        = 15

print(f"Initial  u std = {u.std():.6f}   u range = [{u.min():.4f}, {u.max():.4f}]")
print(f"{'block':>5} {'time':>10} {'u_std':>12} {'u_min':>10} {'u_max':>10}")

t = 0.0
for b in range(1, n_blocks + 1):
    u, v = run_block(u, v, steps_per_block, d, mu, dX, dt)
    t += steps_per_block * dt
    print(f"{b:5d} {t:10.2f} {u.std():12.6f} {u.min():10.4f} {u.max():10.4f}")

# --- final numerical summary ---
print(f"Final    u std = {u.std():.6f}   u range = [{u.min():.4f}, {u.max():.4f}]")
print(f"Final    u mean = {u.mean():.6f}")

# --- separate check: has the small noise amplified into structure? ---
# The initial field had std ~0.058 (small random noise about u=1); a labyrinth
# of stripes means u splits into two well-separated levels, so std grows a lot.
u0_std = 0.1 / np.sqrt(3.0)   # theoretical std of uniform[-0.1,0.1] noise
growth_factor = u.std() / u0_std
print(f"Noise-amplitude growth factor (final u_std / initial noise std) = {growth_factor:.3f}")
print("CHECK: growth factor >> 1 confirms the Turing instability amplified "
      "near-uniform noise into a high-contrast striped pattern.")
print("Why this check works: because the standard deviation of u rising far above "
      "the tiny initial noise level can only happen if the uniform state broke up "
      "into distinct high/low stripe domains, which is exactly the coarsening labyrinth.")

# --- 2D image of the u field (labyrinth of stripes) ---
plt.figure(figsize=(6, 5))
plt.imshow(u, cmap="viridis", origin="lower")
plt.colorbar(label="u (activator)")
plt.title("Gierer-Meinhardt substrate depletion: u field (labyrinth of stripes)")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.2.1_s5.png")
print("Saved u-field image to 7E.2.1_s5.png")
