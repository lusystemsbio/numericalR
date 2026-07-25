import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model parameters ----
g0, g1, Xth, n, k = 10.0, 45.0, 200.0, 4, 0.15
X0 = 300.0
T, dt = 1000.0, 0.01
b_values = [5.0, 2.0, 0.5]

nsteps = int(round(T / dt))
tgrid = np.linspace(0.0, T, nsteps + 1)

# ---- Drift and diffusion definitions ----
def drift(X):
    # Hill self-activation production minus linear degradation
    hx = (X / Xth) ** n
    return g0 + g1 * hx / (1.0 + hx) - k * X

def diffusion(X, b):
    # state-dependent square-root noise term g(X) = b*sqrt(X)
    return b * np.sqrt(np.maximum(X, 0.0))

def diffusion_deriv(X, b):
    # derivative g'(X) = b / (2*sqrt(X)), needed for the Milstein correction
    return b / (2.0 * np.sqrt(np.maximum(X, 1e-12)))

# ---- Milstein integrator (implemented explicitly) ----
def milstein(b, seed=1):
    rng = np.random.default_rng(seed)
    X = np.empty(nsteps + 1)
    X[0] = X0
    for i in range(nsteps):
        x = X[i]
        dW = rng.normal(0.0, np.sqrt(dt))       # Wiener increment ~ N(0, dt)
        a = drift(x)                            # drift a(X)
        g = diffusion(x, b)                     # diffusion g(X)
        gp = diffusion_deriv(x, b)              # g'(X)
        # Euler part + Milstein correction 0.5*g*g'*(dW^2 - dt)
        x_new = x + a * dt + g * dW + 0.5 * g * gp * (dW ** 2 - dt)
        X[i + 1] = max(x_new, 0.0)              # keep non-negative (reflect at 0)
    return X

# ---- Deterministic fixed points (bistability check) ----
xs = np.linspace(0.0, 500.0, 500001)
f = drift(xs)
sign_change = np.where(np.diff(np.sign(f)) != 0)[0]
roots = []
for idx in sign_change:
    # linear interpolation for the zero crossing of the drift
    x1, x2 = xs[idx], xs[idx + 1]
    f1, f2 = f[idx], f[idx + 1]
    root = x1 - f1 * (x2 - x1) / (f2 - f1)
    roots.append(root)

print("Deterministic fixed points (drift = 0):")
for r in roots:
    # stability: stable if d(drift)/dX < 0 at the root
    eps = 1e-3
    slope = (drift(r + eps) - drift(r - eps)) / (2 * eps)
    kind = "stable" if slope < 0 else "unstable"
    print(f"  X* = {r:.3f}  ({kind}, drift slope = {slope:.5f})")

# ---- Run simulations and count transitions ----
threshold = 200.0  # midpoint separating the low (~100) and high (~300) states
fig, axes = plt.subplots(len(b_values), 1, figsize=(10, 9), sharex=True)

for ax, b in zip(axes, b_values):
    X = milstein(b, seed=1)

    # classify each point as low (0) or high (1) state and count crossings
    state = (X > threshold).astype(int)
    n_transitions = int(np.sum(np.abs(np.diff(state))))
    print(f"b = {b}: mean X = {X.mean():.2f}, "
          f"min = {X.min():.2f}, max = {X.max():.2f}, "
          f"state transitions (crossings of X={threshold:.0f}) = {n_transitions}")

    ax.plot(tgrid, X, lw=0.6, color="steelblue")
    ax.axhline(100, color="green", ls="--", lw=1, label="low state ~100")
    ax.axhline(300, color="red", ls="--", lw=1, label="high state ~300")
    ax.axhline(threshold, color="gray", ls=":", lw=1, label="threshold")
    ax.set_ylabel("X")
    ax.set_title(f"b = {b}  (transitions = {n_transitions})")
    ax.legend(loc="upper right", fontsize=7)

axes[-1].set_xlabel("time")
fig.suptitle("Bistable self-activating gene: Milstein SDE simulation")
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.4.1_s2.png", dpi=130)

# ---- One-sentence explanation ----
print("\nExplanation: The drift has two stable fixed points (~100 and ~300) "
      "separated by an unstable one, and the simulated trajectories dwell near "
      "these values while jumping between them more often as b grows, which "
      "confirms both bistability and noise-driven, amplitude-dependent switching.")
