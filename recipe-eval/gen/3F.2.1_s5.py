import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ----------------------------------------------------------------------
# Vector field f(state) = [dX/dt, dY/dt]
# Each gene: basal + repressive Hill of the OTHER gene - linear degradation
# ----------------------------------------------------------------------
def f(state):
    X, Y = state
    dX = gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X
    dY = gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y
    return np.array([dX, dY])

# Wrappers for solve_ivp (t, y) signature; forward and time-reversed flows
def rhs_forward(t, s):
    return f(s)

def rhs_reversed(t, s):
    return -f(s)   # time-reversed system: dX/dt = -f(X)

# ----------------------------------------------------------------------
# 1) Locate the fixed points, in particular the saddle near (190, 127)
# ----------------------------------------------------------------------
def find_fp(guess):
    sol = fsolve(lambda s: f(s), guess, full_output=True)
    root, info, ier, msg = sol
    return np.array(root), ier == 1

# Three fixed points of a bistable toggle: two stable nodes + one saddle
seeds_fp = [(50, 300), (350, 50), (190, 127)]
fixed_points = []
for g in seeds_fp:
    p, ok = find_fp(g)
    if ok:
        # de-duplicate
        if not any(np.allclose(p, q, atol=1e-3) for q in fixed_points):
            fixed_points.append(p)

print("Fixed points found:")
for p in fixed_points:
    print(f"  ({p[0]:.6f}, {p[1]:.6f})")

# The saddle is the one nearest the supplied guess (190, 127)
saddle_guess = np.array([190.0, 127.0])
saddle, ok = find_fp(saddle_guess)
print(f"Saddle point: ({saddle[0]:.6f}, {saddle[1]:.6f})")

# ----------------------------------------------------------------------
# 2) Jacobian at the saddle -> eigen-decomposition
#    The separatrix is the STABLE manifold of the saddle (eigenvalue < 0).
# ----------------------------------------------------------------------
def jacobian(state, h=1e-6):
    # finite-difference Jacobian of f
    J = np.zeros((2, 2))
    for j in range(2):
        dp = np.zeros(2); dp[j] = h
        J[:, j] = (f(state + dp) - f(state - dp)) / (2 * h)
    return J

J = jacobian(saddle)
evals, evecs = np.linalg.eig(J)
print("Jacobian eigenvalues at saddle:")
for i, lam in enumerate(evals):
    print(f"  lambda_{i} = {lam.real:.6f}  (eigvec = [{evecs[0,i].real:.4f}, {evecs[1,i].real:.4f}])")

# Stable eigen-direction = eigenvector of the NEGATIVE eigenvalue.
# In the time-reversed flow this direction becomes UNSTABLE, so points
# seeded along it are carried outward and trace the separatrix.
i_stable = int(np.argmin(evals.real))
v_stable = evecs[:, i_stable].real
v_stable = v_stable / np.linalg.norm(v_stable)
print(f"Stable eigen-direction (traced by reversed flow): "
      f"[{v_stable[0]:.6f}, {v_stable[1]:.6f}]")

# ----------------------------------------------------------------------
# 3) Seed ten points just beside the saddle and integrate the
#    time-reversed system so the flow carries them along the separatrix.
# ----------------------------------------------------------------------
eps = 0.5                       # tiny offset from the saddle
t_span = (0.0, 4000.0)
t_eval = np.linspace(*t_span, 4000)

sep_branches = []
n_seed = 10
for k in range(n_seed):
    # spread seeds symmetrically about the saddle along +/- stable direction
    sign = 1.0 if k % 2 == 0 else -1.0
    mag = eps * (1 + k // 2)     # a few offset magnitudes on each side
    seed = saddle + sign * mag * v_stable
    sol = solve_ivp(rhs_reversed, t_span, seed, t_eval=t_eval,
                    rtol=1e-9, atol=1e-9, max_step=1.0)
    sep_branches.append(sol.y)

# Report endpoints reached by the two separatrix branches
for k in (0, 1):
    xe, ye = sep_branches[k][0, -1], sep_branches[k][1, -1]
    print(f"Separatrix branch seed#{k} endpoint: ({xe:.3f}, {ye:.3f})")

# ----------------------------------------------------------------------
# 4) Sample forward trajectories from a grid of starts (for the plot)
# ----------------------------------------------------------------------
def integrate_forward(start, T=2000.0):
    sol = solve_ivp(rhs_forward, (0, T), start,
                    t_eval=np.linspace(0, T, 2000),
                    rtol=1e-8, atol=1e-8)
    return sol.y

# Identify the two stable attractors (exclude the saddle)
attractors = [p for p in fixed_points
              if not np.allclose(p, saddle, atol=1e-2)]
if len(attractors) < 2:
    # fall back: integrate far-flung starts to recover attractors
    a1 = integrate_forward([20, 400])[:, -1]
    a2 = integrate_forward([400, 20])[:, -1]
    attractors = [a1, a2]
print("Stable attractors:")
for a in attractors:
    print(f"  ({a[0]:.4f}, {a[1]:.4f})")

sample_starts = [(60, 350), (100, 300), (250, 200),
                 (300, 100), (350, 60), (150, 250)]
forward_samples = [integrate_forward(np.array(s, float)) for s in sample_starts]

# ----------------------------------------------------------------------
# 5) CHECK: random forward starts should each end in one basin and
#    never cross the separatrix. We classify by final attractor and
#    verify separation across the traced curve.
# ----------------------------------------------------------------------
rng = np.random.default_rng(0)
n_random = 200
def which_basin(endpoint):
    d = [np.linalg.norm(endpoint - a) for a in attractors]
    return int(np.argmin(d))

# Build a single monotone separatrix polyline (both branches + saddle)
b0 = sep_branches[0]      # one side of stable direction
b1 = sep_branches[1]      # other side
sep_x = np.concatenate([b0[0][::-1], [saddle[0]], b1[0]])
sep_y = np.concatenate([b0[1][::-1], [saddle[1]], b1[1]])
order = np.argsort(sep_x)
sep_x_sorted, sep_y_sorted = sep_x[order], sep_y[order]

# Signed side of a point relative to the separatrix (interpolate Y at its X)
def side_of_separatrix(pt):
    x, y = pt
    if x < sep_x_sorted[0] or x > sep_x_sorted[-1]:
        return None  # outside the traced X-range; skip
    y_sep = np.interp(x, sep_x_sorted, sep_y_sorted)
    return np.sign(y - y_sep)

basin_counts = [0, 0]
side_by_basin = {0: [], 1: []}
starts = rng.uniform([0, 0], [450, 500], size=(n_random, 2))
for s in starts:
    end = integrate_forward(s, T=3000.0)[:, -1]
    b = which_basin(end)
    basin_counts[b] += 1
    sd = side_of_separatrix(s)
    if sd is not None and sd != 0:
        side_by_basin[b].append(sd)

print(f"Random forward starts: basin0={basin_counts[0]}, basin1={basin_counts[1]}")

# Consistency: within the separatrix's X-range, each basin should sit
# predominantly on ONE side of the curve (no crossings).
def purity(sides):
    if not sides:
        return float('nan')
    sides = np.array(sides)
    return max((sides > 0).mean(), (sides < 0).mean())

p0 = purity(side_by_basin[0])
p1 = purity(side_by_basin[1])
print(f"Basin0 one-sidedness fraction: {p0:.4f}")
print(f"Basin1 one-sidedness fraction: {p1:.4f}")
# Check the two basins fall on OPPOSITE sides
mean0 = np.mean(side_by_basin[0]) if side_by_basin[0] else 0.0
mean1 = np.mean(side_by_basin[1]) if side_by_basin[1] else 0.0
print(f"Mean side sign basin0={mean0:.3f}, basin1={mean1:.3f} "
      f"(opposite signs => curve separates basins)")
separates = (p0 > 0.95 and p1 > 0.95 and np.sign(mean0) != np.sign(mean1))
print(f"Separatrix cleanly divides the two basins: {separates}")

# ----------------------------------------------------------------------
# 6) Phase-plane plot
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 7))

# direction field (light)
xx, yy = np.meshgrid(np.linspace(0, 450, 22), np.linspace(0, 500, 22))
U = gX0 + gX1 / (1 + (yy / Yth) ** nY) - kX * xx
V = gY0 + gY1 / (1 + (xx / Xth) ** nX) - kY * yy
N = np.hypot(U, V); N[N == 0] = 1
ax.quiver(xx, yy, U / N, V / N, color="0.8", pivot="mid", scale=40, width=0.002)

# forward sample trajectories
for i, tr in enumerate(forward_samples):
    ax.plot(tr[0], tr[1], color="tab:blue", lw=1.2,
            label="forward trajectory" if i == 0 else None)

# separatrix branches (time-reversed flow)
for i, br in enumerate(sep_branches[:2]):
    ax.plot(br[0], br[1], color="crimson", lw=2.5,
            label="separatrix (stable manifold)" if i == 0 else None)

# fixed points
ax.plot(saddle[0], saddle[1], "ks", ms=9, label="saddle")
for i, a in enumerate(attractors):
    ax.plot(a[0], a[1], "g*", ms=15,
            label="stable attractor" if i == 0 else None)

ax.set_xlim(0, 450); ax.set_ylim(0, 500)
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Toggle-switch separatrix via time-reversed integration")
ax.legend(loc="upper right", fontsize=8)
fig.tight_layout()
fig.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3F.2.1_s5.png", dpi=130)

# ----------------------------------------------------------------------
# Why the check works (one sentence):
print("Explanation: Because the separatrix is the stable manifold of the "
      "saddle, forward trajectories are repelled from it and can only "
      "converge to the attractor on their own side; so if every random "
      "start ends in the basin matching the side it began on and none "
      "cross the traced curve, that curve is exactly the basin boundary.")
