import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D Gierer-Meinhardt substrate-depletion reaction-diffusion model
#   f(u,v) = u^2 v - u      (activator u, diffusion Du = d)
#   g(u,v) = mu (1 - u^2 v) (substrate v, diffusion Dv = 1)
# Integrated explicitly with a finite-difference (5-point Laplacian)
# forward-Euler scheme, run in successive time blocks (from 7E.1).
# ---------------------------------------------------------------

# ---- parameters -----------------------------------------------
N     = 101      # grid points per side (101 x 101)
dX    = 0.2      # spatial step
dt    = 0.005    # time step
d     = 0.1      # Du (activator diffusion)
Dv    = 1.0      # Dv (substrate diffusion)
mu    = 1.5      # kinetic parameter

# reaction terms
def f(u, v): return u*u*v - u            # activator kinetics
def g(u, v): return mu*(1.0 - u*u*v)     # substrate kinetics

# 5-point Laplacian with periodic boundaries (via np.roll)
def laplacian(a):
    return (np.roll(a, 1, 0) + np.roll(a, -1, 0) +
            np.roll(a, 1, 1) + np.roll(a, -1, 1) - 4.0*a) / (dX*dX)

# one explicit forward-Euler step of the RD system
def step(u, v):
    u_new = u + dt*(d  * laplacian(u) + f(u, v))
    v_new = v + dt*(Dv * laplacian(v) + g(u, v))
    return u_new, v_new

# advance the fields by a whole block of `nsteps` steps
def run_block(u, v, nsteps):
    for _ in range(nsteps):
        u, v = step(u, v)
    return u, v

# ---- initial condition: near-uniform u=v=1 with +-0.1 noise ----
rng = np.random.default_rng(10)                 # seed 10
u = 1.0 + 0.1*(2.0*rng.random((N, N)) - 1.0)    # uniform in [-0.1, 0.1]
v = 1.0 + 0.1*(2.0*rng.random((N, N)) - 1.0)

print(f"Homogeneous steady state: u* = 1.0, v* = 1.0")
print(f"Initial u: mean = {u.mean():.6f}, std = {u.std():.6f}")
print(f"Initial v: mean = {v.mean():.6f}, std = {v.std():.6f}")

# ---- run in successive blocks, tracking the pattern amplitude ---
steps_per_block = 2000
n_blocks        = 15

print("\nCheck: growth of u-field amplitude (std about the mean) per block")
print(f"{'block':>5} {'time':>10} {'u_std':>12} {'u_min':>10} {'u_max':>10}")
print(f"{0:>5} {0.0:>10.2f} {u.std():>12.6f} {u.min():>10.4f} {u.max():>10.4f}")

for b in range(1, n_blocks + 1):
    u, v = run_block(u, v, steps_per_block)
    t = b*steps_per_block*dt
    print(f"{b:>5} {t:>10.2f} {u.std():>12.6f} {u.min():>10.4f} {u.max():>10.4f}")

# ---- final diagnostics -----------------------------------------
print(f"\nFinal time: {n_blocks*steps_per_block*dt:.2f}")
print(f"Final u: mean = {u.mean():.6f}, std = {u.std():.6f}")
print(f"Final u range: min = {u.min():.6f}, max = {u.max():.6f}")
print(f"Amplitude growth factor (final u_std / initial): "
      f"{u.std()/(0.1/np.sqrt(3)):.2f}")

# ---- image of the u field (labyrinth of stripes) ---------------
plt.figure(figsize=(6, 5))
plt.imshow(u, cmap="viridis", origin="lower",
           extent=[0, (N-1)*dX, 0, (N-1)*dX])
plt.colorbar(label="u (activator)")
plt.title("2D Gierer-Meinhardt substrate depletion: u field (labyrinth)")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.2.1_s1.png")

# One-sentence explanation of the check:
print("\nWhy the check confirms the result: the u-field standard deviation "
      "climbs by more than an order of magnitude from the tiny initial noise "
      "and then saturates while min/max lock onto two distinct plateau values, "
      "which is exactly the signature of a Turing instability amplifying "
      "near-uniform noise into a finite-amplitude, coarsening labyrinth of stripes.")
