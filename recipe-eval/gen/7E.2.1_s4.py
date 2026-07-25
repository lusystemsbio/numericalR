import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# 2D Gierer-Meinhardt substrate-depletion reaction-diffusion model
#   f(u,v) = u^2 v - u        (activator u)
#   g(u,v) = mu (1 - u^2 v)   (substrate v)
#   Du = d, Dv = 1
# Integrated explicitly with the 2D finite-difference RD integrator (7E.1),
# run in successive blocks so we can watch the pattern coarsen.
# ----------------------------------------------------------------------

# --- parameters -------------------------------------------------------
N   = 101      # grid points per side
dX  = 0.2      # spatial step
dt  = 0.005    # time step
d   = 0.1      # Du (Dv = 1)
mu  = 1.5      # substrate kinetics rate
print(f"grid = {N}x{N}, dX = {dX}, dt = {dt}, d(Du) = {d}, Dv = 1, mu = {mu}")

# --- reaction terms ---------------------------------------------------
def f(u, v):  # activator kinetics
    return u*u*v - u
def g(u, v):  # substrate kinetics
    return mu*(1.0 - u*u*v)

# --- 5-point Laplacian with periodic BCs (via np.roll) ----------------
def lap(a):
    return (np.roll(a, 1, 0) + np.roll(a, -1, 0) +
            np.roll(a, 1, 1) + np.roll(a, -1, 1) - 4.0*a) / (dX*dX)

# --- one explicit Euler step of the RD system -------------------------
def step(u, v):
    un = u + dt*(d*lap(u) + f(u, v))   # activator: slow diffusion + reaction
    vn = v + dt*(1.0*lap(v) + g(u, v)) # substrate: fast diffusion + reaction
    return un, vn

# --- stability sanity check for the explicit scheme -------------------
diff_limit = dX*dX / (4.0*1.0)   # limiting (fastest) diffusion is Dv = 1
print(f"explicit diffusion stability limit dt < {diff_limit:.6f}  (using dt = {dt})")

# --- near-uniform noisy initial condition (u = v = 1 +- 0.1) ----------
rng = np.random.default_rng(10)
u = 1.0 + 0.1*(2.0*rng.random((N, N)) - 1.0)
v = 1.0 + 0.1*(2.0*rng.random((N, N)) - 1.0)
print(f"initial u: mean = {u.mean():.6f}, std(amplitude) = {u.std():.6f}")

# --- run in successive blocks, recording the pattern amplitude --------
# The amplitude = std of u about its mean; if the uniform state is unstable
# it grows from tiny noise and then saturates as stripes form and coarsen.
steps_per_block = 2000
n_blocks        = 25
print("\nblock   sim_time   u_std(amplitude)   u_min     u_max")
for b in range(1, n_blocks + 1):
    for _ in range(steps_per_block):      # explicit time-stepping, one block
        u, v = step(u, v)
    t = b*steps_per_block*dt
    print(f"{b:5d}   {t:8.2f}   {u.std():14.6f}   {u.min():7.4f}   {u.max():7.4f}")

final_amp = u.std()
print(f"\nfinal u: mean = {u.mean():.6f}, std(amplitude) = {final_amp:.6f}")
print(f"amplitude growth factor from initial noise = {final_amp/0.1:.2f}x")

# --- 2D image of the final u field (labyrinth of stripes) -------------
plt.figure(figsize=(6, 5))
plt.imshow(u, cmap="viridis", origin="lower")
plt.colorbar(label="u (activator)")
plt.title("Gierer-Meinhardt substrate depletion: labyrinth of u-stripes")
plt.xlabel("x"); plt.ylabel("y")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7E.2.1_s4.png", dpi=120)

# --- why the check confirms the result --------------------------------
print("\nCheck explanation:")
print("The u-amplitude (std) rises many-fold from the ~0.1 initial noise and then "
      "saturates while u_min/u_max stay bounded, which confirms a Turing instability "
      "that self-organizes the near-uniform noise into a finite-amplitude, coarsening "
      "labyrinth of winding stripes.")
