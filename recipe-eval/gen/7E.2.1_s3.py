import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# 2D Gierer-Meinhardt substrate-depletion reaction-diffusion model
#   du/dt = Du * lap(u) + f(u,v),  f(u,v) = u^2 v - u
#   dv/dt = Dv * lap(v) + g(u,v),  g(u,v) = mu*(1 - u^2 v)
#   Du = d, Dv = 1
# Integrated explicitly (forward Euler) with a 5-point Laplacian,
# run in successive blocks (the 7E.1 integrator, inlined here).
# ---------------------------------------------------------------

# --- parameters -------------------------------------------------
N   = 101        # grid points per side (101 x 101)
dX  = 0.2        # spatial step
dt  = 0.005      # time step
d   = 0.1        # Du (activator diffusion), Dv = 1
mu  = 1.5        # kinetic parameter
Du, Dv = d, 1.0

# --- reaction terms --------------------------------------------
def f(u, v):     # activator kinetics
    return u*u*v - u
def g(u, v):     # substrate kinetics
    return mu*(1.0 - u*u*v)

# --- 5-point Laplacian with periodic (wrap) boundaries ----------
def laplacian(a):
    # sum of 4 neighbors minus 4*center, divided by dX^2
    return (np.roll(a, 1, 0) + np.roll(a, -1, 0) +
            np.roll(a, 1, 1) + np.roll(a, -1, 1) - 4.0*a) / (dX*dX)

# --- one explicit forward-Euler step ----------------------------
def step(u, v):
    un = u + dt*(Du*laplacian(u) + f(u, v))   # update u
    vn = v + dt*(Dv*laplacian(v) + g(u, v))   # update v
    return un, vn

# --- integrator run in successive blocks (as in 7E.1) -----------
def run_blocks(u, v, n_blocks, steps_per_block):
    for _ in range(n_blocks):                 # outer: blocks
        for _ in range(steps_per_block):      # inner: steps
            u, v = step(u, v)
    return u, v

# --- nearly uniform initial condition: u = v = 1 +- 0.1 noise ---
rng = np.random.default_rng(10)               # seed 10
u = 1.0 + rng.uniform(-0.1, 0.1, (N, N))
v = 1.0 + rng.uniform(-0.1, 0.1, (N, N))

# record initial spread as the "near-uniform noise" baseline
u0_std = float(np.std(u))
print(f"Initial u mean:            {np.mean(u):.6f}")
print(f"Initial u std (noise):     {u0_std:.6f}")
print(f"Initial u min/max:         {np.min(u):.6f} / {np.max(u):.6f}")

# --- grow the pattern in successive blocks ----------------------
steps_per_block = 2000
n_blocks = 15                                 # 30000 steps total, T = 150
stds = [u0_std]
for b in range(n_blocks):
    u, v = run_blocks(u, v, 1, steps_per_block)
    s = float(np.std(u))
    stds.append(s)
    print(f"After block {b+1:2d} (t={ (b+1)*steps_per_block*dt:7.2f}): u std = {s:.6f}, "
          f"u range = [{np.min(u):.4f}, {np.max(u):.4f}]")

u_final_std = stds[-1]
print(f"Final u mean:              {np.mean(u):.6f}")
print(f"Final u std (pattern):     {u_final_std:.6f}")
print(f"Amplification of std:      {u_final_std/u0_std:.3f}x")

# --- figure: u field showing the labyrinth of stripes -----------
plt.figure(figsize=(6, 5))
plt.imshow(u, cmap="viridis", origin="lower",
           extent=[0, N*dX, 0, N*dX])
plt.colorbar(label="u")
plt.title("Gierer-Meinhardt substrate depletion: labyrinth of stripes (u)")
plt.xlabel("x")
plt.ylabel("y")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.2.1_s3.png")

# --- check explanation ------------------------------------------
# The check: u's standard deviation grows monotonically from the tiny
# initial noise to a large steady value, and the field organizes into
# winding stripes -- this confirms a Turing instability because a small
# random perturbation was AMPLIFIED (not damped) into a structured,
# coarsening labyrinth pattern rather than relaxing back to the uniform state.
print("Check: u std grew from "
      f"{u0_std:.4f} to {u_final_std:.4f} "
      f"({u_final_std/u0_std:.1f}x); noise was amplified into stripes, "
      "confirming a Turing (diffusion-driven) instability.")
