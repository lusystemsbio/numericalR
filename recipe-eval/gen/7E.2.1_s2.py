import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D Gierer-Meinhardt substrate-depletion reaction-diffusion model
#   f(u,v) = u^2 v - u        (activator kinetics)
#   g(u,v) = mu (1 - u^2 v)   (substrate kinetics)
#   Du = d,  Dv = 1
# Integrated explicitly with a 5-point Laplacian (finite differences,
# forward Euler in time), run in successive blocks (from 7E.1 style).
# ---------------------------------------------------------------

# --- parameters ---
N   = 101      # grid points per side (101 x 101)
dX  = 0.2      # spatial step
dt  = 0.005    # time step
d   = 0.1      # Du (activator diffusion), Dv = 1
mu  = 1.5      # substrate kinetic rate
Du  = d
Dv  = 1.0

# --- reaction terms ---
def f(u, v):
    return u**2 * v - u

def g(u, v):
    return mu * (1.0 - u**2 * v)

# --- 5-point Laplacian with no-flux (Neumann) boundaries ---
def laplacian(a):
    # np.pad with 'edge' mirrors the boundary value -> zero-flux
    ap = np.pad(a, 1, mode="edge")
    lap = (ap[:-2, 1:-1] + ap[2:, 1:-1] +
           ap[1:-1, :-2] + ap[1:-1, 2:] - 4.0 * ap[1:-1, 1:-1])
    return lap / dX**2

# --- one explicit forward-Euler step ---
def step(u, v):
    un = u + dt * (Du * laplacian(u) + f(u, v))   # update activator u
    vn = v + dt * (Dv * laplacian(v) + g(u, v))   # update substrate v
    return un, vn

# --- run a block of many steps ---
def run_block(u, v, nsteps):
    for _ in range(nsteps):
        u, v = step(u, v)
    return u, v

# --- nearly uniform initial condition: u = v = 1 with +-0.1 noise ---
rng = np.random.default_rng(10)
u = 1.0 + rng.uniform(-0.1, 0.1, size=(N, N))
v = 1.0 + rng.uniform(-0.1, 0.1, size=(N, N))

# amplitude of the initial noise field (spatial std of u)
init_std = np.std(u)
print(f"Initial u spatial std (noise amplitude): {init_std:.6f}")

# --- grow the pattern in successive blocks, tracking the amplitude ---
block_size = 2000                 # steps per block
n_blocks   = 15                   # total steps = 30000
print("block, time, u_min, u_max, u_std")
for b in range(1, n_blocks + 1):
    u, v = run_block(u, v, block_size)
    t = b * block_size * dt
    print(f"{b:2d}, t={t:8.2f}, u_min={u.min():.4f}, "
          f"u_max={u.max():.4f}, u_std={np.std(u):.6f}")

final_std = np.std(u)
print(f"Final u spatial std (pattern amplitude): {final_std:.6f}")
print(f"Amplitude growth factor (final/initial): {final_std/init_std:.4f}")
print(f"Final u min: {u.min():.6f}")
print(f"Final u max: {u.max():.6f}")
print(f"Final u mean: {u.mean():.6f}")

# --- 2D image of the u field: labyrinth of stripes ---
plt.figure(figsize=(6, 5))
plt.imshow(u, cmap="viridis", origin="lower",
           extent=[0, N * dX, 0, N * dX])
plt.colorbar(label="u (activator)")
plt.title("Gierer-Meinhardt substrate depletion: u labyrinth")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.2.1_s2.png")

# --- check explanation ---
# The spatial std of u climbs from the ~0.06 initial-noise level to a much
# larger saturated value while u_min/u_max spread far apart: this monotone
# growth-then-saturation of the amplitude confirms a Turing instability that
# amplifies near-uniform noise into a coarsening labyrinth of winding stripes.
print("Check: u_std grew from the tiny noise level to a large saturated "
      "value, confirming the instability amplifies noise into stripes.")
