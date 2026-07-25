import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Gierer-Meinhardt activator-inhibitor Turing pattern (1D)
# f(u,v) = u^2/v - u   (activator: local self-enhancement)
# g(u,v) = mu*(u^2 - v) (inhibitor: long-range, fast diffusion)
# Du = d (slow), Dv = 1 (fast)
# ---------------------------------------------------------------

# --- parameters ---
L    = 20.0     # domain length
dX   = 0.2      # spatial step
dt   = 0.01     # time step
d    = 0.1      # activator diffusion (slow)
Dv   = 1.0      # inhibitor diffusion (fast)
mu   = 1.5      # inhibitor kinetic rate

N  = int(round(L / dX)) + 1     # number of grid points
X  = np.linspace(0.0, L, N)

# --- reaction kinetics ---
def f(u, v):
    return u**2 / v - u          # activator kinetics

def g(u, v):
    return mu * (u**2 - v)       # inhibitor kinetics

# --- 1D Laplacian with zero-flux (Neumann) boundaries ---
def laplacian(c):
    lap = np.empty_like(c)
    # interior points: standard central difference
    lap[1:-1] = (c[2:] - 2.0*c[1:-1] + c[:-2]) / dX**2
    # boundaries: reflect (ghost node = neighbour) -> no flux
    lap[0]  = (2.0*c[1]  - 2.0*c[0])  / dX**2
    lap[-1] = (2.0*c[-2] - 2.0*c[-1]) / dX**2
    return lap

# --- explicit one-step forward-Euler update of the RD system ---
def step(u, v):
    u_new = u + dt * (d  * laplacian(u) + f(u, v))
    v_new = v + dt * (Dv * laplacian(v) + g(u, v))
    return u_new, v_new

# --- integrate a block of many steps (the multi-component integrator, 7C.1) ---
def run_block(u, v, n_steps):
    for _ in range(n_steps):
        u, v = step(u, v)
    return u, v

# --- nearly uniform initial condition u = v = 1 with +-0.1 noise ---
rng = np.random.default_rng(10)
u = 1.0 + 0.1 * (2.0 * rng.random(N) - 1.0)   # uniform noise in [-0.1, 0.1]
v = 1.0 + 0.1 * (2.0 * rng.random(N) - 1.0)

# --- stability (explicit diffusion number must be < 0.5) ---
diff_number = Dv * dt / dX**2
print(f"Explicit diffusion number (Dv*dt/dX^2): {diff_number:.4f}")

# --- run in successive time blocks until the pattern is stationary ---
steps_per_block = 20000            # 200 time units per block
n_blocks        = 15               # total 3000 time units
prev_u = u.copy()
for b in range(n_blocks):
    u, v = run_block(u, v, steps_per_block)
    block_change = np.max(np.abs(u - prev_u))
    print(f"Block {b+1:2d}  time={ (b+1)*steps_per_block*dt:7.1f}  max|u change| since last block: {block_change:.3e}")
    prev_u = u.copy()

# ---------------------------------------------------------------
# CHECK: stationary periodic Turing pattern with u,v IN PHASE
# ---------------------------------------------------------------

# 1) Stationarity: run one more block and measure how little u moves.
u2, v2 = run_block(u, v, steps_per_block)
stationarity = np.max(np.abs(u2 - u))
print(f"\nStationarity residual (max|u| change over extra block): {stationarity:.3e}")

# 2) Periodicity: count interior peaks of u to confirm a periodic pattern.
def count_peaks(c):
    return int(np.sum((c[1:-1] > c[:-2]) & (c[1:-1] > c[2:])))
n_peaks = count_peaks(u)
print(f"Number of interior peaks in u (periodicity): {n_peaks}")

# 3) Phase relation: spatial correlation between u and v.
#    Positive correlation  => u and v peak together (activator-inhibitor).
#    Negative correlation  => out of phase (substrate depletion).
uc = u - np.mean(u)
vc = v - np.mean(v)
corr = np.sum(uc * vc) / np.sqrt(np.sum(uc**2) * np.sum(vc**2))
print(f"Spatial correlation between u and v: {corr:.4f}")

# 4) Confirm peak locations coincide.
u_peak_idx = np.argmax(u)
v_peak_idx = np.argmax(v)
print(f"Location of global u maximum: X = {X[u_peak_idx]:.2f}")
print(f"Location of global v maximum: X = {X[v_peak_idx]:.2f}")
print(f"u ranges [{u.min():.3f}, {u.max():.3f}], v ranges [{v.min():.3f}, {v.max():.3f}]")

in_phase = corr > 0
print(f"Pattern is IN PHASE (u and v peak together): {in_phase}")
# One-sentence explanation of the check:
print("Check meaning: a POSITIVE spatial u-v correlation with coincident maxima and a"
      " negligible stationarity residual confirms a fixed periodic Turing pattern whose"
      " activator and inhibitor peak together, unlike the anti-phase substrate-depletion case.")

# --- plot the stationary pattern ---
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(X, u, 'b-o', ms=3, label='activator u(X)')
ax.plot(X, v, 'r-s', ms=3, label='inhibitor v(X)')
ax.set_xlabel('X')
ax.set_ylabel('concentration')
ax.set_title('Gierer-Meinhardt Turing pattern: u and v peak together (in phase)')
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7C.3.1_s5.png")
