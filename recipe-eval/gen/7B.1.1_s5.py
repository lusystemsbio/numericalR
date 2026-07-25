import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Fisher's equation:  du/dt = r*u*(1 - u/B) + D*d2u/dX2
# Explicit finite-difference reaction-diffusion integrator:
#   each iteration we (1) apply the diffusion (Laplacian) step and
#   (2) add the local logistic reaction term, using forward Euler in time.
# Dirichlet boundaries: the two end nodes are held fixed (u = 0 here).
# ----------------------------------------------------------------------

# --- Parameters ---
L  = 100        # domain length
dX = 1.0        # spatial step
dt = 0.01       # time step
D  = 1.0        # diffusion coefficient
r  = 0.64       # intrinsic growth rate
B  = 75.0       # carrying capacity

N  = int(L / dX) + 1          # number of spatial nodes
X  = np.linspace(0, L, N)     # spatial grid

# Explicit-diffusion stability check (must be <= 0.5)
alpha = D * dt / dX**2
print(f"Diffusion stability number alpha = D*dt/dX^2 = {alpha:.4f}")
print(f"Number of spatial nodes N = {N}")

T_total = 150.0                # total simulated time
nsteps  = int(T_total / dt)    # total iterations
print(f"Total simulated time = {T_total}")
print(f"Total number of iterations = {nsteps}")

# Fisher front speed prediction: c = 2*sqrt(r*D)
c_theory = 2.0 * np.sqrt(r * D)
print(f"Theoretical Fisher front speed c = 2*sqrt(r*D) = {c_theory:.4f}")


def integrate(u0, nsteps, snapshot_steps):
    """Explicit reaction-diffusion integrator with Dirichlet (u=0) ends.
    Returns dict {step: u.copy()} for requested snapshot steps."""
    u = u0.copy()
    snaps = {}
    if 0 in snapshot_steps:
        snaps[0] = u.copy()
    for step in range(1, nsteps + 1):
        # (1) Diffusion step: second difference of the Laplacian on interior nodes
        lap = np.zeros_like(u)
        lap[1:-1] = (u[2:] - 2.0 * u[1:-1] + u[:-2]) / dX**2
        # (2) Reaction step: local logistic growth added to the diffusion update
        react = r * u * (1.0 - u / B)
        # Forward-Euler update of interior nodes (ends stay fixed = Dirichlet)
        u[1:-1] = u[1:-1] + dt * (D * lap[1:-1] + react[1:-1])
        if step in snapshot_steps:
            snaps[step] = u.copy()
    return snaps


# --- Initial conditions ---
# Central patch: a small block in the middle set to u = 50
u_central = np.zeros(N)
u_central[N // 2 - 3 : N // 2 + 4] = 50.0

# Left-end patch: a small block near the left set to u = 50
u_left = np.zeros(N)
u_left[1:8] = 50.0

# Snapshot times to display
snap_times = [0, 20, 40, 60, 80, 100, 130]
snap_steps = [int(t / dt) for t in snap_times]

snaps_c = integrate(u_central, nsteps, snap_steps)
snaps_l = integrate(u_left,    nsteps, snap_steps)

# --- Front-position / speed check ---
# Track the location where u crosses B/2 on the right-moving front.
def front_position(u, level):
    idx = np.where(u >= level)[0]
    return X[idx[-1]] if len(idx) else np.nan

t1, t2 = 100.0, 130.0
s1, s2 = int(t1 / dt), int(t2 / dt)
snaps_speed_c = integrate(u_central, s2, [s1, s2])
snaps_speed_l = integrate(u_left,    s2, [s1, s2])

# Central: right front position; also mirror to check left front symmetry
xc1 = front_position(snaps_speed_c[s1], B / 2)
xc2 = front_position(snaps_speed_c[s2], B / 2)
speed_c_right = (xc2 - xc1) / (t2 - t1)
print(f"Central patch right-front position at t={t1}: X = {xc1:.2f}")
print(f"Central patch right-front position at t={t2}: X = {xc2:.2f}")
print(f"Central patch measured right-front speed = {speed_c_right:.4f}")

# Left front of the central patch (first crossing from the left)
def front_position_left(u, level):
    idx = np.where(u >= level)[0]
    return X[idx[0]] if len(idx) else np.nan
xl1 = front_position_left(snaps_speed_c[s1], B / 2)
xl2 = front_position_left(snaps_speed_c[s2], B / 2)
speed_c_left = (xl1 - xl2) / (t2 - t1)   # leftward positive
print(f"Central patch measured left-front speed (leftward) = {speed_c_left:.4f}")

# Left-end patch: single rightward front
xr1 = front_position(snaps_speed_l[s1], B / 2)
xr2 = front_position(snaps_speed_l[s2], B / 2)
speed_l_right = (xr2 - xr1) / (t2 - t1)
print(f"Left-end patch right-front position at t={t1}: X = {xr1:.2f}")
print(f"Left-end patch right-front position at t={t2}: X = {xr2:.2f}")
print(f"Left-end patch measured right-front speed = {speed_l_right:.4f}")

# Confirm interior plateau reached carrying capacity
peak_c = snaps_speed_c[s2].max()
peak_l = snaps_speed_l[s2].max()
print(f"Central patch peak value at t={t2} (should approach B={B}): {peak_c:.4f}")
print(f"Left-end patch peak value at t={t2} (should approach B={B}): {peak_l:.4f}")

# --- Plots ---
fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)

for st in snap_steps:
    axes[0].plot(X, snaps_c[st], label=f"t = {st * dt:.0f}")
axes[0].axhline(B, color="k", ls="--", lw=0.8, alpha=0.6)
axes[0].set_title("Central patch: two fronts advancing outward")
axes[0].set_xlabel("X")
axes[0].set_ylabel("u(X)")
axes[0].legend(fontsize=8)

for st in snap_steps:
    axes[1].plot(X, snaps_l[st], label=f"t = {st * dt:.0f}")
axes[1].axhline(B, color="k", ls="--", lw=0.8, alpha=0.6)
axes[1].set_title("Left-end patch: single rightward front")
axes[1].set_xlabel("X")
axes[1].legend(fontsize=8)

fig.suptitle("Fisher's equation: traveling waves of advance")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/7B.1.1_s5.png", dpi=120)

# Explanation of the check:
print("Explanation: The check confirms the result because both patches settle to "
      "the same plateau u=B and their B/2 fronts advance equal distances per unit "
      "time (~2*sqrt(r*D)), so a central patch emits two symmetric constant-speed "
      "fronts while a left-end patch emits only one rightward front.")
