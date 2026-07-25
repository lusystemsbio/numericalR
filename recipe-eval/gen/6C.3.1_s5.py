import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(12345)

# ----------------------------------------------------------------------
# Milstein integrator for the scalar SDE
#     dX = f(X) dt + s(X) dW,     s(X) = sqrt(2 D(X))
# implemented explicitly (no black-box SDE solver).
#
# Euler-Maruyama step : X + f*dt + s*dW
# Milstein correction : + 0.5 * s * s' * (dW^2 - dt)
#   -> this stochastic-Taylor term vanishes when s'(X) == 0 (constant noise).
# ----------------------------------------------------------------------
def milstein_trajectory(f, s, sprime, x0, T, N, dW=None, rng=rng):
    dt = T / N                       # fixed time step
    if dW is None:                   # Wiener increments dW ~ N(0, dt)
        dW = rng.normal(0.0, np.sqrt(dt), size=N)
    X = np.empty(N + 1)
    X[0] = x0
    for n in range(N):
        x = X[n]
        drift = f(x) * dt            # deterministic part
        diff = s(x) * dW[n]          # Euler-Maruyama stochastic part
        corr = 0.5 * s(x) * sprime(x) * (dW[n] ** 2 - dt)  # Milstein correction
        X[n + 1] = x + drift + diff + corr
    return X, dW

def euler_maruyama_trajectory(f, s, x0, T, N, dW=None, rng=rng):
    # Same as Milstein but with the correction term dropped.
    dt = T / N
    if dW is None:
        dW = rng.normal(0.0, np.sqrt(dt), size=N)
    X = np.empty(N + 1)
    X[0] = x0
    for n in range(N):
        x = X[n]
        X[n + 1] = x + f(x) * dt + s(x) * dW[n]
    return X, dW

# ----------------------------------------------------------------------
# (1) Integrator applied to a test SDE -> returns a trajectory.
#     Gene-circuit model of 6C.4: self-activating gene with birth-death
#     (state-dependent, multiplicative) noise.
#         f(X) = a + b * X^h/(K^h + X^h) - g*X     (Hill activation - decay)
#         D(X) = X   ->  s(X) = sqrt(2X),  s'(X) = 1/sqrt(2X)
# ----------------------------------------------------------------------
a, b, K, h, g = 0.5, 4.0, 2.0, 4.0, 1.0

def f_gene(x):
    xp = max(x, 0.0)
    return a + b * xp**h / (K**h + xp**h) - g * xp

def s_gene(x):
    return np.sqrt(2.0 * max(x, 0.0))

def sprime_gene(x):
    # d/dx sqrt(2x) = 1/sqrt(2x)
    xp = max(x, 1e-12)
    return 1.0 / np.sqrt(2.0 * xp)

T_gene, N_gene, x0_gene = 20.0, 4000, 1.0
traj_gene, _ = milstein_trajectory(f_gene, s_gene, sprime_gene, x0_gene, T_gene, N_gene)
t_gene = np.linspace(0.0, T_gene, N_gene + 1)

print("Gene-circuit trajectory (Milstein):")
print("  initial X       =", traj_gene[0])
print("  final X         =", traj_gene[-1])
print("  trajectory min  =", traj_gene.min())
print("  trajectory max  =", traj_gene.max())
print("  trajectory mean =", traj_gene.mean())

# ----------------------------------------------------------------------
# (2) Convergence check.
#     Use a model with a KNOWN exact solution so we can measure strong error.
#
#  (a) STATE-DEPENDENT noise: geometric Brownian motion
#         dX = mu X dt + sig X dW,   s'(X) = sig   (nonzero, X-dependent)
#      exact: X_T = X0 * exp((mu - sig^2/2) T + sig W_T)
#      Expectation: Milstein strong order ~1.0, Euler-Maruyama ~0.5.
#
#  (b) CONSTANT noise: dX = mu X dt + sig dW,   s'(X) = 0
#      Milstein correction term is identically zero, so Milstein == EM.
# ----------------------------------------------------------------------
def strong_convergence_slopes(f, s, sprime, exact_fn, x0, T, Ns, n_paths, rng):
    em_err, mil_err = [], []
    for N in Ns:
        dt = T / N
        e_em = e_mil = 0.0
        for _ in range(n_paths):
            dW = rng.normal(0.0, np.sqrt(dt), size=N)  # shared path for both schemes
            W_T = dW.sum()
            x_exact = exact_fn(x0, T, W_T)
            x_em, _ = euler_maruyama_trajectory(f, s, x0, T, N, dW=dW)
            x_mil, _ = milstein_trajectory(f, s, sprime, x0, T, N, dW=dW)
            e_em += abs(x_em[-1] - x_exact)
            e_mil += abs(x_mil[-1] - x_exact)
        em_err.append(e_em / n_paths)
        mil_err.append(e_mil / n_paths)
    dts = np.array([T / N for N in Ns])
    slope_em = np.polyfit(np.log(dts), np.log(em_err), 1)[0]
    slope_mil = np.polyfit(np.log(dts), np.log(mil_err), 1)[0]
    return dts, np.array(em_err), np.array(mil_err), slope_em, slope_mil

Ns = [16, 32, 64, 128, 256]
n_paths = 400
T_c, x0_c, mu, sig = 1.0, 1.0, 1.0, 0.5

# (a) state-dependent (multiplicative) noise
f_gbm = lambda x: mu * x
s_gbm = lambda x: sig * x
sprime_gbm = lambda x: sig
exact_gbm = lambda x0, T, W: x0 * np.exp((mu - 0.5 * sig**2) * T + sig * W)

dts, em_a, mil_a, sl_em_a, sl_mil_a = strong_convergence_slopes(
    f_gbm, s_gbm, sprime_gbm, exact_gbm, x0_c, T_c, Ns, n_paths, np.random.default_rng(1))

# (b) constant noise
f_add = lambda x: mu * x
s_add = lambda x: sig
sprime_add = lambda x: 0.0
def exact_add(x0, T, W):
    # dX = mu X dt + sig dW  ->  X_T = x0 e^{mu T} + sig * integral, but the
    # discretized reference below uses the same shared-path Ito integral; for
    # slope purposes we use a fine EM reference is not needed since s'=0 makes
    # Milstein and EM byte-identical. We compare the two schemes directly.
    return None

# For constant noise we directly verify the two schemes coincide.
rng_b = np.random.default_rng(2)
max_diff = 0.0
for N in Ns:
    dt = T_c / N
    dW = rng_b.normal(0.0, np.sqrt(dt), size=N)
    x_em, _ = euler_maruyama_trajectory(f_add, s_add, x0_c, T_c, N, dW=dW)
    x_mil, _ = milstein_trajectory(f_add, s_add, sprime_add, x0_c, T_c, N, dW=dW)
    max_diff = max(max_diff, np.max(np.abs(x_em - x_mil)))

print()
print("Strong-convergence order (state-dependent GBM noise):")
print("  Euler-Maruyama slope =", sl_em_a, "(expected ~0.5)")
print("  Milstein slope       =", sl_mil_a, "(expected ~1.0)")
print("  EM   errors          =", em_a)
print("  Mil  errors          =", mil_a)
print()
print("Constant-noise check (correction term should vanish):")
print("  max |Milstein - EulerMaruyama| over all steps/grids =", max_diff)
print()
print("Why this confirms the result:")
print("  The Milstein correction lowers strong error to O(dt) (slope ~1) versus")
print("  O(sqrt(dt)) (slope ~0.5) for Euler-Maruyama ONLY when s'(X)!=0; for")
print("  constant noise s'(X)=0 makes the correction exactly zero, so the two")
print("  integrators produce byte-identical trajectories (max diff ~= 0),")
print("  proving the improvement is due solely to state-dependent noise.")

# ----------------------------------------------------------------------
# Plot: gene-circuit trajectory + convergence comparison.
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(t_gene, traj_gene, lw=0.8, color="C0")
ax1.set_xlabel("time")
ax1.set_ylabel("X")
ax1.set_title("Milstein trajectory: gene-circuit SDE (6C.4)")

ax2.loglog(dts, em_a, "o-", label=f"Euler-Maruyama (slope {sl_em_a:.2f})")
ax2.loglog(dts, mil_a, "s-", label=f"Milstein (slope {sl_mil_a:.2f})")
ax2.loglog(dts, dts**0.5 * em_a[0] / dts[0]**0.5, "k--", lw=0.8, label="ref slope 0.5")
ax2.loglog(dts, dts**1.0 * mil_a[0] / dts[0]**1.0, "k:", lw=0.8, label="ref slope 1.0")
ax2.set_xlabel("dt")
ax2.set_ylabel("mean strong error at T")
ax2.set_title("Strong convergence (state-dependent noise)")
ax2.legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.3.1_s5.png")
print()
print("Saved figure to 6C.3.1_s5.png")
