import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(12345)

# ----------------------------------------------------------------------
# One explicit step of each integrator for dX = f(X) dt + s(X) dW.
# We write the stochastic-Taylor terms out by hand (no black-box routine).
# ----------------------------------------------------------------------
def euler_maruyama_step(x, f, s, sprime, dt, dW):
    # Euler-Maruyama: deterministic drift + first-order noise term.
    return x + f(x) * dt + s(x) * dW

def milstein_step(x, f, s, sprime, dt, dW):
    # Milstein = Euler-Maruyama + stochastic-Taylor correction.
    em = x + f(x) * dt + s(x) * dW            # base Euler-Maruyama step
    correction = 0.5 * s(x) * sprime(x) * (dW**2 - dt)  # vanishes when s'(x)=0
    return em + correction

# ----------------------------------------------------------------------
# Integrate an SDE, returning the full trajectory. If dW_path is given we
# reuse those Wiener increments (dW ~ N(0, dt)); otherwise we draw them.
# ----------------------------------------------------------------------
def integrate_sde(step_fn, f, s, sprime, x0, dt, n_steps, dW_path=None, gen=rng):
    if dW_path is None:
        dW_path = gen.normal(0.0, np.sqrt(dt), size=n_steps)  # dW ~ N(0, dt)
    traj = np.empty(n_steps + 1)
    traj[0] = x0
    x = x0
    for i in range(n_steps):
        x = step_fn(x, f, s, sprime, dt, dW_path[i])
        traj[i + 1] = x
    return traj, dW_path

# ----------------------------------------------------------------------
# (1) Apply the integrator to the gene-circuit model of 6C.4.
#     Chemical-Langevin form: production - degradation, with
#     state-dependent noise s(X) = sqrt(2 D(X)), D(X) = 0.5*(k + gamma*X).
# ----------------------------------------------------------------------
k, gamma = 5.0, 1.0
f_gene      = lambda x: k - gamma * x
D_gene      = lambda x: 0.5 * (k + gamma * max(x, 0.0))
s_gene      = lambda x: np.sqrt(2.0 * D_gene(x))
sprime_gene = lambda x: gamma / (2.0 * np.sqrt(2.0 * D_gene(x)))  # d/dx sqrt(k+gamma*x)

dt_gene, T_gene = 0.001, 20.0
n_gene = int(T_gene / dt_gene)
gene_traj, _ = integrate_sde(milstein_step, f_gene, s_gene, sprime_gene,
                             x0=k / gamma, dt=dt_gene, n_steps=n_gene)

print("Gene-circuit trajectory (Milstein):")
print(f"  number of steps            = {n_gene}")
print(f"  final state X(T)           = {gene_traj[-1]:.6f}")
print(f"  trajectory mean            = {gene_traj.mean():.6f}")
print(f"  trajectory std            = {gene_traj.std():.6f}")
print(f"  deterministic steady state = {k / gamma:.6f}")

# ----------------------------------------------------------------------
# (2) Strong-convergence study on geometric Brownian motion (state-dependent
#     noise, exact solution known). We reuse ONE fine Wiener path per sample
#     and coarsen it, comparing each integrator's endpoint to the exact one.
# ----------------------------------------------------------------------
mu, sigma, x0_gbm, T = 0.5, 0.4, 1.0, 1.0
f_gbm      = lambda x: mu * x
s_gbm      = lambda x: sigma * x
sprime_gbm = lambda x: sigma           # s'(x) = sigma  (noise depends on state)

n_paths = 400
dt_fine = T / 2048
n_fine = int(T / dt_fine)
levels = [2**e for e in range(1, 7)]   # coarsening factors -> dt = level*dt_fine

err_em = np.zeros(len(levels))
err_mil = np.zeros(len(levels))
for _ in range(n_paths):
    dW_fine = rng.normal(0.0, np.sqrt(dt_fine), size=n_fine)
    W_T = dW_fine.sum()
    x_exact = x0_gbm * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)  # exact endpoint
    for j, L in enumerate(levels):
        dt = L * dt_fine
        n_steps = n_fine // L
        # aggregate fine increments into coarse ones (same Brownian path)
        dW = dW_fine[: n_steps * L].reshape(n_steps, L).sum(axis=1)
        em_traj, _ = integrate_sde(euler_maruyama_step, f_gbm, s_gbm, sprime_gbm,
                                   x0_gbm, dt, n_steps, dW_path=dW)
        mil_traj, _ = integrate_sde(milstein_step, f_gbm, s_gbm, sprime_gbm,
                                    x0_gbm, dt, n_steps, dW_path=dW)
        err_em[j] += abs(em_traj[-1] - x_exact)
        err_mil[j] += abs(mil_traj[-1] - x_exact)
err_em /= n_paths
err_mil /= n_paths

dts = np.array(levels) * dt_fine
slope_em = np.polyfit(np.log(dts), np.log(err_em), 1)[0]
slope_mil = np.polyfit(np.log(dts), np.log(err_mil), 1)[0]

print("\nStrong convergence on GBM (state-dependent noise):")
for dt, ee, em in zip(dts, err_em, err_mil):
    print(f"  dt={dt:.5f}  EM error={ee:.6e}  Milstein error={em:.6e}")
print(f"  Euler-Maruyama strong order (slope) = {slope_em:.4f}")
print(f"  Milstein strong order (slope)       = {slope_mil:.4f}")

# ----------------------------------------------------------------------
# (3) Constant-noise check: additive SDE, s'(x)=0, so the correction is 0
#     and Milstein must reproduce Euler-Maruyama exactly.
# ----------------------------------------------------------------------
f_add      = lambda x: mu * x
s_add      = lambda x: sigma          # constant noise
sprime_add = lambda x: 0.0            # s'(x) = 0
dt_c, n_c = 0.01, 100
dW_c = rng.normal(0.0, np.sqrt(dt_c), size=n_c)
em_c, _ = integrate_sde(euler_maruyama_step, f_add, s_add, sprime_add, 1.0, dt_c, n_c, dW_path=dW_c)
mil_c, _ = integrate_sde(milstein_step, f_add, s_add, sprime_add, 1.0, dt_c, n_c, dW_path=dW_c)
max_diff = np.max(np.abs(em_c - mil_c))

print("\nConstant-noise check (additive SDE, s'(x)=0):")
print(f"  max |Milstein - EM| over trajectory = {max_diff:.3e}")
print("  Explanation: because the correction 0.5*s*s'*(dW^2-dt) is identically "
      "zero when s'(x)=0, Milstein and Euler-Maruyama coincide for constant noise, "
      "so the correction can only help when the noise depends on the state (as shown "
      "by Milstein's order ~1 vs EM's order ~0.5 for GBM).")

# ----------------------------------------------------------------------
# Figure: convergence (state-dependent) + gene trajectory.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
ax[0].loglog(dts, err_em, "o-", label=f"Euler-Maruyama (slope {slope_em:.2f})")
ax[0].loglog(dts, err_mil, "s-", label=f"Milstein (slope {slope_mil:.2f})")
ax[0].loglog(dts, err_mil[0] * (dts / dts[0]), "k--", alpha=0.5, label="order 1 ref")
ax[0].loglog(dts, err_em[0] * (dts / dts[0])**0.5, "k:", alpha=0.5, label="order 1/2 ref")
ax[0].set_xlabel("dt"); ax[0].set_ylabel("strong error at T")
ax[0].set_title("State-dependent noise (GBM)"); ax[0].legend(fontsize=8)

t_gene = np.linspace(0, T_gene, n_gene + 1)
ax[1].plot(t_gene, gene_traj, lw=0.8)
ax[1].axhline(k / gamma, color="r", ls="--", label="det. steady state")
ax[1].set_xlabel("t"); ax[1].set_ylabel("X(t)")
ax[1].set_title("Gene-circuit trajectory (Milstein)"); ax[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.3.1_s1.png")
