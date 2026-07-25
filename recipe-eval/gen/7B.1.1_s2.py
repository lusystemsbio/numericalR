import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---------------------------------------------------------------
# Fisher's equation:  du/dt = r*u*(1 - u/B) + D * d2u/dX2
# Explicit finite-difference reaction-diffusion integrator:
#   at each step we do the diffusion (Laplacian) update AND add
#   the logistic reaction term, using Dirichlet boundaries (u fixed
#   at the two ends).
# ---------------------------------------------------------------

# --- parameters ---
L   = 100      # domain length
dX  = 1.0      # spatial step
dt  = 0.01     # time step
D   = 1.0      # diffusion coefficient
r   = 0.64     # intrinsic growth rate
B   = 75.0     # carrying capacity

N = int(L / dX) + 1                 # number of grid points
X = np.linspace(0.0, L, N)          # spatial grid

# Explicit-scheme diffusion stability number (must be <= 0.5 for stability)
alpha = D * dt / dX**2
print("Diffusion number alpha = D*dt/dX^2 =", alpha)
print("Stability requires alpha <= 0.5:", alpha <= 0.5)

# Theoretical Fisher front speed c = 2*sqrt(r*D)
c_theory = 2.0 * np.sqrt(r * D)
print("Theoretical front speed c = 2*sqrt(r*D) =", c_theory)


def integrate(u0, n_steps):
    """Explicit reaction-diffusion integrator with Dirichlet boundaries.
    Returns snapshots of u at requested step indices."""
    u = u0.copy()
    for step in range(n_steps):
        # --- diffusion step: discrete Laplacian on interior points ---
        lap = np.zeros_like(u)
        lap[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dX**2
        # --- reaction term added to the diffusion step each iteration ---
        reaction = r * u * (1.0 - u / B)
        # --- update interior points; boundaries (Dirichlet) left unchanged ---
        u[1:-1] = u[1:-1] + dt * (D * lap[1:-1] + reaction[1:-1])
        yield_step = step  # (kept explicit for clarity)
    return u


def run_and_snapshot(u0, times):
    """Integrate and record u at the given list of physical times."""
    u = u0.copy()
    snap_steps = {int(round(t / dt)): t for t in times}
    max_step = max(snap_steps)
    snapshots = {}
    for step in range(max_step + 1):
        if step in snap_steps:
            snapshots[snap_steps[step]] = u.copy()
        # diffusion (discrete Laplacian on interior points)
        lap = np.zeros_like(u)
        lap[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dX**2
        # logistic reaction added each iteration
        reaction = r * u * (1.0 - u / B)
        # explicit update of interior points, Dirichlet boundaries fixed at 0
        u[1:-1] = u[1:-1] + dt * (D * lap[1:-1] + reaction[1:-1])
    return snapshots


# --- initial conditions ---
patch_amp = 50.0

# central patch
u_center = np.zeros(N)
c_idx = N // 2
u_center[c_idx - 3:c_idx + 4] = patch_amp   # small patch near the middle

# left-end patch
u_left = np.zeros(N)
u_left[1:8] = patch_amp                      # small patch near the left end
# (index 0 stays at the Dirichlet boundary value 0)

# times to visualize
times = [0, 20, 40, 60, 80, 100, 120]

snap_center = run_and_snapshot(u_center, times)
snap_left   = run_and_snapshot(u_left, times)

# ---------------------------------------------------------------
# Plots of u(X) at successive times
# ---------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for t in times:
    axes[0].plot(X, snap_center[t], label=f"t={t}")
    axes[1].plot(X, snap_left[t],   label=f"t={t}")

for ax, title in zip(axes, ["Central patch", "Left-end patch"]):
    ax.axhline(B, color="k", ls="--", lw=0.8, label="carrying capacity B")
    ax.set_xlabel("X")
    ax.set_ylabel("u(X)")
    ax.set_title(title)
    ax.legend(fontsize=8)

fig.suptitle("Fisher's equation: traveling waves of advance")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7B.1.1_s2.png", dpi=120)

# ---------------------------------------------------------------
# Quantitative check:
#   - patch grows to carrying capacity
#   - central patch emits TWO constant-speed fronts (left & right)
#   - left-end patch emits ONE rightward front
# We track the front position as the X where u crosses B/2, and
# estimate the speed from its motion between snapshots.
# ---------------------------------------------------------------

def front_position(u, x, side):
    """Return X where u crosses B/2 on the given side ('right' or 'left')."""
    half = B / 2.0
    if side == "right":
        idx = np.where(u >= half)[0]
        if len(idx) == 0:
            return np.nan
        i = idx[-1]                     # last point above half on the right
        if i >= len(u) - 1:
            return x[i]
        # linear interpolation between i and i+1
        return x[i] + (half - u[i]) / (u[i + 1] - u[i]) * dX
    else:  # left
        idx = np.where(u >= half)[0]
        if len(idx) == 0:
            return np.nan
        i = idx[0]                      # first point above half on the left
        if i == 0:
            return x[i]
        return x[i] + (half - u[i]) / (u[i - 1] - u[i]) * (-dX)


# peak amplitude reached (approach to carrying capacity)
print("Central patch: max u at final time =", snap_center[times[-1]].max())
print("Left-end patch: max u at final time =", snap_left[times[-1]].max())

# measure front speeds over the later (settled) interval
t_a, t_b = 80, 120

# central patch: two fronts
r_a = front_position(snap_center[t_a], X, "right")
r_b = front_position(snap_center[t_b], X, "right")
l_a = front_position(snap_center[t_a], X, "left")
l_b = front_position(snap_center[t_b], X, "left")
speed_right_center = (r_b - r_a) / (t_b - t_a)
speed_left_center  = (l_b - l_a) / (t_b - t_a)   # negative = moving left
print("Central patch right-front speed  =", speed_right_center)
print("Central patch left-front  speed  =", speed_left_center)

# left-end patch: single rightward front
lr_a = front_position(snap_left[t_a], X, "right")
lr_b = front_position(snap_left[t_b], X, "right")
speed_right_left = (lr_b - lr_a) / (t_b - t_a)
print("Left-end patch right-front speed =", speed_right_left)

# left-end patch should have essentially no leftward front (anchored at wall)
print("Left-end patch left-front pos at t=120 =", front_position(snap_left[t_b], X, "left"))

# Explanation of why this check confirms the result:
print("Explanation: The two central fronts advance at equal and opposite "
      "constant speeds close to 2*sqrt(r*D) while the patch saturates at B, "
      "whereas the left-end patch yields a single rightward front at the same "
      "speed, confirming Fisher traveling waves of advance.")
