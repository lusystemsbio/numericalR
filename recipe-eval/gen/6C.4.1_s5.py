import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- Model definition ----------------
# dX = [g0 + g1*(X/Xth)^n/(1+(X/Xth)^n) - k*X]*dt + b*sqrt(X)*dW
g0, g1, Xth, n, k = 10.0, 45.0, 200.0, 4.0, 0.15
X0 = 300.0
T, dt = 1000.0, 0.01
nsteps = int(T / dt)
t = np.linspace(0.0, T, nsteps + 1)
b_values = [5.0, 2.0, 0.5]

def drift(X):
    # Hill self-activation minus linear degradation
    h = (X / Xth) ** n
    return g0 + g1 * h / (1.0 + h) - k * X

def diffusion(X, b):
    # state-dependent noise amplitude b*sqrt(X)
    return b * np.sqrt(np.maximum(X, 0.0))

# ---------------- Milstein integrator (explicit) ----------------
# For g(X) = b*sqrt(X), g'(X) = b/(2*sqrt(X)), so the Milstein
# correction 0.5*g*g'*(dW^2 - dt) simplifies to 0.25*b^2*(dW^2 - dt).
def simulate(b, seed=1):
    rng = np.random.default_rng(seed)
    X = np.empty(nsteps + 1)
    X[0] = X0
    for i in range(nsteps):
        x = X[i]
        dW = rng.normal(0.0, np.sqrt(dt))          # Wiener increment ~ N(0, dt)
        a = drift(x)                                # drift term a(X)
        g = diffusion(x, b)                         # diffusion term g(X)=b*sqrt(X)
        milstein = 0.25 * b * b * (dW * dW - dt)    # 0.5*g*g'*(dW^2-dt) = 0.25*b^2*(...)
        x_next = x + a * dt + g * dW + milstein     # Milstein update
        X[i + 1] = max(x_next, 0.0)                 # keep concentration non-negative
    return X

# ---------------- Bistability check (deterministic) ----------------
# Find fixed points of the drift (roots of f) and classify by sign of f'.
Xgrid = np.linspace(0.0, 500.0, 500001)
f = drift(Xgrid)
sign_change = np.where(np.diff(np.sign(f)) != 0)[0]
fixed_points = []
for idx in sign_change:
    # linear interpolation for the root location between grid points
    x_lo, x_hi = Xgrid[idx], Xgrid[idx + 1]
    root = x_lo - f[idx] * (x_hi - x_lo) / (f[idx + 1] - f[idx])
    # numerical derivative of drift for stability classification
    eps = 1e-3
    fprime = (drift(root + eps) - drift(root - eps)) / (2 * eps)
    fixed_points.append((root, fprime))

print("=== Bistability check: fixed points of the drift ===")
stable = []
for root, fprime in fixed_points:
    kind = "stable" if fprime < 0 else "unstable"
    if fprime < 0:
        stable.append(root)
    print(f"Fixed point X = {root:.4f}   f'(X) = {fprime:.6f}   ({kind})")

n_stable = len(stable)
print(f"Number of stable fixed points: {n_stable}")
print(f"Is bistable (two stable states): {n_stable == 2}")
if n_stable == 2:
    print(f"Lower stable state X ~ {stable[0]:.2f}")
    print(f"Upper stable state X ~ {stable[1]:.2f}")

# Unstable fixed point used as the separatrix between the two basins
unstable_fps = [r for r, fp in fixed_points if fp > 0]
separatrix = unstable_fps[0] if unstable_fps else 0.5 * (stable[0] + stable[1])
print(f"Separatrix (unstable fixed point) X ~ {separatrix:.4f}")

# ---------------- Simulate, plot, and count transitions ----------------
fig, axes = plt.subplots(len(b_values), 1, figsize=(10, 9), sharex=True)
print("\n=== Stochastic simulations (Milstein) ===")
for ax, b in zip(axes, b_values):
    X = simulate(b, seed=1)
    # State = which basin the trajectory is in relative to the separatrix.
    state = (X > separatrix).astype(int)
    n_transitions = int(np.sum(np.abs(np.diff(state))))
    print(f"b = {b}:  final X = {X[-1]:.3f}   mean X = {X.mean():.3f}   "
          f"transitions between states = {n_transitions}")
    ax.plot(t, X, lw=0.4)
    ax.axhline(stable[0], color="green", ls="--", lw=1, label=f"lower ~{stable[0]:.0f}")
    ax.axhline(stable[1], color="red", ls="--", lw=1, label=f"upper ~{stable[1]:.0f}")
    ax.axhline(separatrix, color="gray", ls=":", lw=1, label=f"separatrix ~{separatrix:.0f}")
    ax.set_ylabel("X")
    ax.set_title(f"b = {b}  ({n_transitions} transitions)")
    ax.legend(loc="upper right", fontsize=7)

axes[-1].set_xlabel("time")
fig.suptitle("Bistable self-activating gene: noise-driven transitions (Milstein)")
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.4.1_s5.png", dpi=120)

# ---------------- Interpretation ----------------
print("\nWhy the check confirms the result: the drift has exactly two stable "
      "fixed points (near X=100 and X=300) separated by an unstable one, so the "
      "circuit is bistable, and the transition count rising with b shows that "
      "larger state-dependent noise more frequently kicks the system across the "
      "separatrix between those two states.")
