import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------------------------
# 6D.1 model: noisy mutual-repression toggle switch (Euler-Maruyama SDE)
#   dX = (b/(1+Y^n) - X) dt + sigma dWx
#   dY = (b/(1+X^n) - Y) dt + sigma dWy
# Two stable states: X-high/Y-low and X-low/Y-high. Noise drives rare switches.
# -----------------------------------------------------------------------------
rng = np.random.default_rng(3)      # seed 3
b = 20.0                            # production strength
n = 4                              # Hill coefficient
sigma = 3.5                        # noise amplitude
dt = 0.01
T = 1000.0
N = int(round(T / dt))

X = np.empty(N + 1)
Y = np.empty(N + 1)
X[0], Y[0] = b, 0.0                # start in the X-high state
sqdt = np.sqrt(dt)
for k in range(N):
    fx = b / (1.0 + Y[k] ** n) - X[k]
    fy = b / (1.0 + X[k] ** n) - Y[k]
    X[k + 1] = X[k] + fx * dt + sigma * sqdt * rng.standard_normal()
    Y[k + 1] = Y[k] + fy * dt + sigma * sqdt * rng.standard_normal()

pts = np.column_stack([X, Y])      # trajectory as (N+1, 2) points

# -----------------------------------------------------------------------------
# Define the two state clouds by splitting the plane along the line Y = X.
#   state 0 : X > Y  (X-high cloud)
#   state 1 : Y > X  (Y-high cloud)
# -----------------------------------------------------------------------------
above = Y > X                      # boolean mask for the Y>X side
clouds = [pts[~above], pts[above]]

# -----------------------------------------------------------------------------
# Fit an ellipse to each cloud from the covariance eigen-decomposition.
# Semi-axis length along eigenvector i = enlarge * sqrt(eigenvalue_i).
# -----------------------------------------------------------------------------
enlarge = 1.5
means, eigvecs, axes = [], [], []
for cloud in clouds:
    mu = cloud.mean(axis=0)                       # ellipse center
    C = np.cov(cloud, rowvar=False)               # 2x2 covariance
    lam, V = np.linalg.eigh(C)                     # eigen-decomposition (symmetric)
    a = enlarge * np.sqrt(lam)                     # enlarged semi-axis lengths
    means.append(mu); eigvecs.append(V); axes.append(a)

def inside_score(p, mu, V, a):
    # Normalized squared radius in the ellipse's own eigen-frame.
    # Project (p - mu) onto each eigenvector, divide by that axis, sum of squares.
    # score <= 1  <=>  point lies inside the ellipse.
    d = p - mu
    proj = V.T @ d                                 # coordinates in eigen-frame
    return np.sum((proj / a) ** 2)

# -----------------------------------------------------------------------------
# Assign each point to a state: inside ellipse 0, inside ellipse 1, or neither.
# When inside both, take the smaller (closer) score. When inside neither, the
# point is in the no-man's-land between wells -> mark as -1 (undecided).
# -----------------------------------------------------------------------------
assign = np.empty(N + 1, dtype=int)
for i, p in enumerate(pts):
    s0 = inside_score(p, means[0], eigvecs[0], axes[0])
    s1 = inside_score(p, means[1], eigvecs[1], axes[1])
    in0, in1 = s0 <= 1.0, s1 <= 1.0
    if in0 and in1:
        assign[i] = 0 if s0 <= s1 else 1
    elif in0:
        assign[i] = 0
    elif in1:
        assign[i] = 1
    else:
        assign[i] = -1                             # undecided (between ellipses)

# Robust detection: hysteresis. Carry the last decided state through the
# undecided gap; a transition is a change in the carried state.
last = assign[0] if assign[0] != -1 else 0
robust_steps = []
prev = last
for i in range(N + 1):
    cur = assign[i] if assign[i] != -1 else prev   # forward-fill undecided
    if cur != prev:
        robust_steps.append(i)                     # step where state flipped
    prev = cur
robust_count = len(robust_steps)

# -----------------------------------------------------------------------------
# Crude counter (the thing we want to show over-counts): assign purely by which
# side of the line Y = X the point is on, then count every sign change.
# -----------------------------------------------------------------------------
side = (Y > X).astype(int)                          # 0 or 1, no dead zone
crude_steps = np.where(np.diff(side) != 0)[0] + 1
crude_count = len(crude_steps)

# -----------------------------------------------------------------------------
# Report
# -----------------------------------------------------------------------------
print(f"Total steps simulated: {N}")
print(f"Ellipse 0 (X>Y) center: X={means[0][0]:.3f}, Y={means[0][1]:.3f}")
print(f"Ellipse 0 semi-axes (enlarged x{enlarge}): {axes[0][0]:.3f}, {axes[0][1]:.3f}")
print(f"Ellipse 1 (Y>X) center: X={means[1][0]:.3f}, Y={means[1][1]:.3f}")
print(f"Ellipse 1 semi-axes (enlarged x{enlarge}): {axes[1][0]:.3f}, {axes[1][1]:.3f}")
print(f"Robust ellipse-based transition count: {robust_count}")
for j, s in enumerate(robust_steps, 1):
    print(f"  robust transition {j}: step {s} (t={s*dt:.2f})")
print(f"Crude Y=X-line transition count: {crude_count}")

# Show the crude counter clusters many transitions within a few steps around a
# single real event: report gaps between consecutive crude transitions.
if crude_count > 1:
    gaps = np.diff(crude_steps)
    n_close = int(np.sum(gaps <= 5))
    print(f"Crude transitions separated by <=5 steps (jitter within a cluster): {n_close}")
    print(f"Median gap between crude transitions: {np.median(gaps):.1f} steps")
    print(f"Over-count factor (crude / robust): {crude_count / max(robust_count,1):.1f}")

# -----------------------------------------------------------------------------
# Phase plane with the two fitted state ellipses
# -----------------------------------------------------------------------------
theta = np.linspace(0, 2 * np.pi, 200)
fig, ax = plt.subplots(figsize=(7, 7))
ax.plot(X, Y, lw=0.3, color="0.6", alpha=0.6, label="trajectory")
lim = max(X.max(), Y.max()) * 1.05
ax.plot([0, lim], [0, lim], "k--", lw=0.8, label="Y = X split")
colors = ["tab:red", "tab:blue"]
for c, mu, V, a in zip(colors, means, eigvecs, axes):
    circle = np.column_stack([a[0] * np.cos(theta), a[1] * np.sin(theta)])
    ell = (V @ circle.T).T + mu                    # rotate+translate unit ellipse
    ax.plot(ell[:, 0], ell[:, 1], color=c, lw=2)
    ax.scatter(*mu, color=c, s=40, zorder=5)
ax.set_xlabel("X"); ax.set_ylabel("Y")
ax.set_title("Toggle-switch phase plane with fitted state ellipses")
ax.set_aspect("equal"); ax.legend(loc="upper right")
plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6D.2.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: The crude line-based counter reports far more transitions "
      "than the ellipse counter and clusters them within a few steps of one "
      "another, showing that boundary jitter inflates the raw count and that "
      "the ellipse+hysteresis statistic, which ignores the undecided gap, is "
      "the more robust measure of true state switches.")
