import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

np.random.seed(0)

# ---------------------------------------------------------------------------
# Milstein integrator for  dX = f(X) dt + s(X) dW,  with  s(X) = sqrt(2 D(X))
# ---------------------------------------------------------------------------
# Implemented explicitly (not via a black-box routine) so every term is visible.
# The Milstein update is the Euler-Maruyama step PLUS the stochastic-Taylor
# correction  0.5 * s(X) * s'(X) * (dW^2 - dt), which vanishes when s'(X)=0.
def milstein_integrate(f, s, sprime, X0, T, dt, dW=None):
    n = int(round(T / dt))                 # number of steps
    X = np.empty(n + 1)
    X[0] = X0
    if dW is None:                         # Wiener increments dW ~ N(0, dt)
        dW = np.sqrt(dt) * np.random.randn(n)
    for i in range(n):
        x = X[i]
        drift = f(x) * dt                  # deterministic Euler part
        diffusion = s(x) * dW[i]           # Euler-Maruyama stochastic part
        correction = 0.5 * s(x) * sprime(x) * (dW[i]**2 - dt)  # Milstein term
        X[i + 1] = x + drift + diffusion + correction
    return X

# Same loop but WITHOUT the correction term = plain Euler-Maruyama (for comparison).
def euler_maruyama_integrate(f, s, X0, T, dt, dW=None):
    n = int(round(T / dt))
    X = np.empty(n + 1)
    X[0] = X0
    if dW is None:
        dW = np.sqrt(dt) * np.random.randn(n)
    for i in range(n):
        x = X[i]
        X[i + 1] = x + f(x) * dt + s(x) * dW[i]
    return X

# ---------------------------------------------------------------------------
# (1) Integrator applied to a test SDE -> returns a trajectory (as in 6C.4)
#     Gene-circuit-style model: bistable drift with state-dependent noise.
# ---------------------------------------------------------------------------
# Drift: production (Hill activation) minus linear degradation.
def gene_f(x):
    return x**2 / (1.0 + x**2) - 0.5 * x
# State-dependent diffusion coefficient D(X) and noise s(X)=sqrt(2 D(X)).
def gene_D(x):
    return 0.05 * (0.1 + x)            # noise grows with abundance
def gene_s(x):
    return np.sqrt(2.0 * gene_D(x))
def gene_sprime(x):                    # d/dx sqrt(2*0.05*(0.1+x)) = 0.05/s(x)
    return 0.05 / gene_s(x)

T, dt = 20.0, 0.005
traj = milstein_integrate(gene_f, gene_s, gene_sprime, X0=0.2, T=T, dt=dt)
t = np.linspace(0.0, T, traj.size)
print("Gene-circuit test trajectory: n_points =", traj.size)
print("Gene-circuit test trajectory: X(0)     =", traj[0])
print("Gene-circuit test trajectory: X(T)     =", traj[-1])
print("Gene-circuit test trajectory: mean X   =", traj.mean())

# ---------------------------------------------------------------------------
# (2) Convergence check: correction helps ONLY for state-dependent noise.
#     Test SDE = geometric Brownian motion, which has an exact solution.
#       dX = mu X dt + sigma X dW,  s(X)=sigma X (state-dependent), s'(X)=sigma
#     Strong error is measured against the exact path built on the SAME dW.
# ---------------------------------------------------------------------------
mu, sigma, X0, Tc = 1.0, 0.8, 1.0, 1.0

def strong_errors(state_dependent):
    if state_dependent:
        f = lambda x: mu * x
        s = lambda x: sigma * x
        sp = lambda x: sigma
    else:  # constant-noise variant: dX = mu X dt + sigma dW, so s'(X)=0
        f = lambda x: mu * x
        s = lambda x: sigma + 0.0 * x
        sp = lambda x: 0.0 * x
    dts = [2.0**-k for k in range(5, 11)]
    M = 400                                  # number of sample paths
    em_err, mil_err = [], []
    for dt in dts:
        n = int(round(Tc / dt))
        e_em = e_mil = 0.0
        for _ in range(M):
            dW = np.sqrt(dt) * np.random.randn(n)
            W_T = dW.sum()                   # Brownian value at T for exact soln
            xem = euler_maruyama_integrate(f, s, X0, Tc, dt, dW)[-1]
            xmil = milstein_integrate(f, s, sp, X0, Tc, dt, dW)[-1]
            if state_dependent:              # exact GBM endpoint on this path
                x_exact = X0 * np.exp((mu - 0.5 * sigma**2) * Tc + sigma * W_T)
            else:                            # exact linear-noise endpoint (Ito)
                x_exact = X0 * np.exp(mu * Tc) + sigma * (
                    np.exp(mu * Tc) * (dW * np.exp(-mu * dt * np.arange(1, n + 1))).sum())
            e_em += abs(xem - x_exact)
            e_mil += abs(xmil - x_exact)
        em_err.append(e_em / M)
        mil_err.append(e_mil / M)
    return np.array(dts), np.array(em_err), np.array(mil_err)

# State-dependent noise: Milstein (order ~1.0) should beat EM (order ~0.5).
dts, em_sd, mil_sd = strong_errors(state_dependent=True)
p_em_sd = np.polyfit(np.log(dts), np.log(em_sd), 1)[0]
p_mil_sd = np.polyfit(np.log(dts), np.log(mil_sd), 1)[0]
print("State-dependent noise: Euler-Maruyama strong order ~", round(p_em_sd, 3))
print("State-dependent noise: Milstein       strong order ~", round(p_mil_sd, 3))
print("State-dependent noise: mean EM error at smallest dt  =", em_sd[-1])
print("State-dependent noise: mean Milstein error at smallest dt =", mil_sd[-1])

# Constant noise: correction term is zero, so the two methods are identical.
f_c = lambda x: -x
s_c = lambda x: sigma + 0.0 * x
sp_c = lambda x: 0.0 * x
dt_c = 0.01
dW_c = np.sqrt(dt_c) * np.random.randn(int(round(Tc / dt_c)))
xem_c = euler_maruyama_integrate(f_c, s_c, X0, Tc, dt_c, dW_c)
xmil_c = milstein_integrate(f_c, s_c, sp_c, X0, Tc, dt_c, dW_c)
max_diff = np.max(np.abs(xem_c - xmil_c))
print("Constant noise: max |Milstein - EM| over path =", max_diff)
print("Constant noise: correction reduces to Euler-Maruyama =", max_diff == 0.0)

# One-sentence explanation of why this check confirms the result:
print("Explanation: because the Milstein correction is proportional to s'(X),",
      "it identically vanishes for constant noise (giving byte-identical EM paths)",
      "yet raises the strong convergence order from ~0.5 to ~1.0 when s'(X) != 0,",
      "confirming the term matters only for state-dependent noise.")

# ---------------------------------------------------------------------------
# Figure: sample trajectory and the convergence comparison.
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
ax[0].plot(t, traj, lw=0.8)
ax[0].set_xlabel("t"); ax[0].set_ylabel("X"); ax[0].set_title("Milstein trajectory (gene circuit)")
ax[1].loglog(dts, em_sd, "o-", label="Euler-Maruyama")
ax[1].loglog(dts, mil_sd, "s-", label="Milstein")
ax[1].set_xlabel("dt"); ax[1].set_ylabel("strong error")
ax[1].set_title("State-dependent noise convergence"); ax[1].legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.3.1_s2.png")
