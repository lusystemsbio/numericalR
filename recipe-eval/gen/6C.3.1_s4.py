import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(12345)

# ----------------------------------------------------------------------
# Milstein integrator (implemented explicitly, step by step).
#
# SDE:  dX = f(X) dt + s(X) dW ,  with  s(X) = sqrt(2 D(X))
#
# Euler-Maruyama step:  X_{n+1} = X_n + f(X)*dt + s(X)*dW
# Milstein adds the stochastic-Taylor correction from Ito's lemma applied
# to the diffusion coefficient:   + 0.5 * s(X) * s'(X) * (dW^2 - dt)
# This extra term captures the leading state-dependence of the noise and
# vanishes identically when s'(X) = 0 (constant noise -> pure EM).
# ----------------------------------------------------------------------
def milstein(f, s, sprime, x0, dt, n_steps, dW=None, use_correction=True):
    """Integrate dX = f(X)dt + s(X)dW and return the trajectory (length n_steps+1)."""
    x = np.empty(n_steps + 1)
    x[0] = x0
    if dW is None:                                   # Wiener increments dW ~ N(0, dt)
        dW = rng.normal(0.0, np.sqrt(dt), size=n_steps)
    for n in range(n_steps):
        xn = x[n]
        drift = f(xn) * dt                           # deterministic (drift) part
        diffusion = s(xn) * dW[n]                     # Euler-Maruyama noise part
        step = xn + drift + diffusion
        if use_correction:                           # Milstein correction term
            step += 0.5 * s(xn) * sprime(xn) * (dW[n] ** 2 - dt)
        x[n + 1] = step
    return x, dW


# ----------------------------------------------------------------------
# Gene-circuit model of 6C.4:  a self-activating gene with birth-death
# (state-dependent, multiplicative-like) noise.
#   f(X)   = a * X^h/(K^h + X^h) + b - g*X      (Hill activation + basal - decay)
#   D(X)   = 0.5 * X       ->     s(X) = sqrt(2 D(X)) = sqrt(X)
#   s'(X)  = 0.5 / sqrt(X)
# ----------------------------------------------------------------------
a, K, h, b, g = 4.0, 2.0, 4.0, 0.2, 1.0

def f_gene(x):
    return a * x**h / (K**h + x**h) + b - g * x

def s_gene(x):
    x = np.maximum(x, 1e-12)                         # keep noise real near 0
    return np.sqrt(x)                                # = sqrt(2 * 0.5 * x)

def sprime_gene(x):
    x = np.maximum(x, 1e-12)
    return 0.5 / np.sqrt(x)

# --- Produce the integrator applied to the test SDE: return a trajectory ---
dt = 0.005
T = 20.0
n_steps = int(T / dt)
t = np.linspace(0.0, T, n_steps + 1)
traj, _ = milstein(f_gene, s_gene, sprime_gene, x0=0.5, dt=dt, n_steps=n_steps)

print("=== Gene-circuit trajectory (6C.4) via Milstein ===")
print("dt                          :", dt)
print("n_steps                     :", n_steps)
print("X(0)                        :", traj[0])
print("X(T) final                  :", traj[-1])
print("trajectory mean             :", np.mean(traj))
print("trajectory std              :", np.std(traj))
print("trajectory min              :", np.min(traj))
print("trajectory max              :", np.max(traj))


# ----------------------------------------------------------------------
# Convergence check: does the correction term help?
# Strong error at final time versus a fine reference path, on the SAME
# Brownian motion, over a range of step sizes.  For state-dependent noise
# Milstein is strong order ~1 while Euler-Maruyama is only order ~0.5.
# For constant noise the correction term is exactly zero, so the two
# methods coincide and both are order ~1.
# ----------------------------------------------------------------------
def strong_error(f, s, sprime, x0, T, dt_list, dt_ref, n_paths=200):
    """Mean |X_method(T) - X_ref(T)| over many Brownian paths, per dt."""
    n_ref = int(round(T / dt_ref))
    err_mil = np.zeros(len(dt_list))
    err_em = np.zeros(len(dt_list))
    for _ in range(n_paths):
        dW_ref = rng.normal(0.0, np.sqrt(dt_ref), size=n_ref)   # one fine path
        x_ref, _ = milstein(f, s, sprime, x0, dt_ref, n_ref, dW=dW_ref, use_correction=True)
        ref = x_ref[-1]
        for i, dt_c in enumerate(dt_list):
            m = int(round(dt_c / dt_ref))                       # coarse:fine ratio
            n_c = int(round(T / dt_c))
            dW_c = dW_ref[:n_c * m].reshape(n_c, m).sum(axis=1)  # coarsen increments
            x_mil, _ = milstein(f, s, sprime, x0, dt_c, n_c, dW=dW_c, use_correction=True)
            x_em, _ = milstein(f, s, sprime, x0, dt_c, n_c, dW=dW_c, use_correction=False)
            err_mil[i] += abs(x_mil[-1] - ref)
            err_em[i] += abs(x_em[-1] - ref)
    return err_mil / n_paths, err_em / n_paths

dt_ref = 2**-13
dt_list = np.array([2**-k for k in range(4, 9)], dtype=float)   # coarse step sizes

# (A) state-dependent noise (gene model)
err_mil_v, err_em_v = strong_error(f_gene, s_gene, sprime_gene, 0.5, 4.0, dt_list, dt_ref)

# (B) constant noise: s(X) = const, s'(X) = 0  -> correction vanishes
sig = 0.7
err_mil_c, err_em_c = strong_error(f_gene, lambda x: sig + 0*x, lambda x: 0*x,
                                   0.5, 4.0, dt_list, dt_ref)

def slope(dt_list, err):
    return np.polyfit(np.log(dt_list), np.log(err), 1)[0]

print("\n=== Convergence check (strong error vs dt) ===")
print("step sizes dt               :", dt_list.tolist())
print("[state-dep] Milstein errors :", err_mil_v.tolist())
print("[state-dep] Euler-M errors  :", err_em_v.tolist())
print("[state-dep] Milstein slope  :", slope(dt_list, err_mil_v))
print("[state-dep] Euler-M slope   :", slope(dt_list, err_em_v))
print("[const-noise] Milstein errs :", err_mil_c.tolist())
print("[const-noise] Euler-M errs  :", err_em_c.tolist())
print("[const-noise] Milstein slope:", slope(dt_list, err_mil_c))
print("[const-noise] Euler-M slope :", slope(dt_list, err_em_c))
print("[const-noise] max |Mil-EM|  :", np.max(np.abs(err_mil_c - err_em_c)),
      "(should be ~0: correction term is identically zero)")

# ----------------------------------------------------------------------
# Why this check confirms the result:
# When the noise is state-dependent the Milstein error decays like dt^1
# while Euler-Maruyama only decays like dt^0.5, whereas for constant noise
# the two methods give identical errors (both order 1) because the
# 0.5*s*s'*(dW^2-dt) correction is exactly zero -- so the improvement is
# demonstrably attributable solely to state-dependent noise.
# ----------------------------------------------------------------------

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(t, traj, lw=0.8)
ax[0].set_xlabel("t"); ax[0].set_ylabel("X(t)")
ax[0].set_title("Gene-circuit SDE trajectory (Milstein)")

ax[1].loglog(dt_list, err_mil_v, "o-", label="Milstein, state-dep noise")
ax[1].loglog(dt_list, err_em_v, "s-", label="Euler-M, state-dep noise")
ax[1].loglog(dt_list, err_mil_c, "^--", label="Milstein, const noise")
ax[1].loglog(dt_list, err_em_c, "v--", label="Euler-M, const noise")
ax[1].loglog(dt_list, dt_list, "k:", label="slope 1 ref")
ax[1].loglog(dt_list, np.sqrt(dt_list) * err_em_v[-1] / np.sqrt(dt_list[-1]),
             "k-.", label="slope 1/2 ref")
ax[1].set_xlabel("dt"); ax[1].set_ylabel("strong error at T")
ax[1].set_title("Strong convergence"); ax[1].legend(fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.3.1_s4.png")
print("\nSaved figure to 6C.3.1_s4.png")
