import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(0)

# ----------------------------------------------------------------------
# Milstein integrator for a 1-D SDE:  dX = f(X) dt + s(X) dW,  dW ~ N(0,dt)
# s(X) = sqrt(2*D(X)) is the state-dependent noise amplitude.
# Milstein adds the stochastic-Taylor correction 0.5*s*s'*(dW^2 - dt)
# to the Euler-Maruyama step; this term is written out explicitly below.
# ----------------------------------------------------------------------
def milstein(f, s, sprime, x0, T, N, dW=None):
    dt = T / N                              # fixed time step
    x = np.empty(N + 1)
    x[0] = x0
    if dW is None:
        dW = rng.normal(0.0, np.sqrt(dt), size=N)   # Wiener increments ~ N(0, dt)
    for n in range(N):
        xn = x[n]
        drift = f(xn) * dt                          # deterministic Euler part
        diffusion = s(xn) * dW[n]                    # Euler-Maruyama noise part
        # Milstein correction from the stochastic Taylor expansion;
        # vanishes identically when s' == 0 (constant noise).
        correction = 0.5 * s(xn) * sprime(xn) * (dW[n] ** 2 - dt)
        x[n + 1] = xn + drift + diffusion + correction
    return dt, x

# Euler-Maruyama = Milstein without the correction term (for comparison).
def euler_maruyama(f, s, x0, T, N, dW=None):
    dt = T / N
    x = np.empty(N + 1)
    x[0] = x0
    if dW is None:
        dW = rng.normal(0.0, np.sqrt(dt), size=N)
    for n in range(N):
        x[n + 1] = x[n] + f(x[n]) * dt + s(x[n]) * dW[n]
    return dt, x

# ----------------------------------------------------------------------
# (1) Apply the integrator to the gene-circuit test SDE of 6C.4.
# Self-activating gene: Hill production + basal - degradation,
# with state-dependent (multiplicative-type) noise D(X) ~ X so more
# molecules -> larger fluctuations. s(X) = sqrt(2*D(X)).
# ----------------------------------------------------------------------
a, K, hill, basal, gamma, Dscale = 4.0, 1.0, 4.0, 0.2, 1.0, 0.05
def f_gene(x):
    xc = max(x, 0.0)
    return a * xc**hill / (K**hill + xc**hill) + basal - gamma * xc
def D_gene(x):
    return Dscale * max(x, 0.0) + 1e-6           # keep D > 0
def s_gene(x):
    return np.sqrt(2.0 * D_gene(x))
def sprime_gene(x):                              # d/dx sqrt(2*Dscale*x + eps)
    return Dscale / np.sqrt(2.0 * D_gene(x))

T_gene, N_gene = 20.0, 4000
dt_g, traj = milstein(f_gene, s_gene, sprime_gene, x0=0.5, T=T_gene, N=N_gene)
t_gene = np.linspace(0.0, T_gene, N_gene + 1)
print("Gene-circuit trajectory: dt =", dt_g)
print("Gene-circuit trajectory: X(0) =", traj[0])
print("Gene-circuit trajectory: X(T) =", traj[-1])
print("Gene-circuit trajectory: mean over t =", np.mean(traj))
print("Gene-circuit trajectory: min, max =", np.min(traj), np.max(traj))

# ----------------------------------------------------------------------
# (2) Convergence check with a known exact solution.
# STATE-DEPENDENT noise: geometric Brownian motion dX = mu X dt + sigma X dW,
#   s(X)=sigma X, s'(X)=sigma, exact X_T = X0*exp((mu-sigma^2/2)T + sigma W_T).
# Strong error is measured against the exact solution driven by the SAME path.
# ----------------------------------------------------------------------
mu, sigma, X0, T = 1.0, 0.8, 1.0, 1.0
f_gbm = lambda x: mu * x
s_gbm = lambda x: sigma * x
sp_gbm = lambda x: sigma
Ns = [16, 32, 64, 128, 256, 512]
paths = 400

def strong_errors(Ns, paths):
    err_mil, err_em = [], []
    for N in Ns:
        dt = T / N
        e_m = e_e = 0.0
        for _ in range(paths):
            dW = rng.normal(0.0, np.sqrt(dt), size=N)
            W_T = dW.sum()
            X_exact = X0 * np.exp((mu - 0.5 * sigma**2) * T + sigma * W_T)
            _, xm = milstein(f_gbm, s_gbm, sp_gbm, X0, T, N, dW=dW)
            _, xe = euler_maruyama(f_gbm, s_gbm, X0, T, N, dW=dW)
            e_m += abs(xm[-1] - X_exact)
            e_e += abs(xe[-1] - X_exact)
        err_mil.append(e_m / paths)
        err_em.append(e_e / paths)
    return np.array(err_mil), np.array(err_em)

err_mil, err_em = strong_errors(Ns, paths)
dts = np.array([T / N for N in Ns])
order_mil = np.polyfit(np.log(dts), np.log(err_mil), 1)[0]
order_em = np.polyfit(np.log(dts), np.log(err_em), 1)[0]
print("State-dependent noise (GBM): dt values =", dts.tolist())
print("State-dependent noise (GBM): Milstein strong errors =", err_mil.tolist())
print("State-dependent noise (GBM): Euler-Maruyama strong errors =", err_em.tolist())
print("State-dependent noise (GBM): Milstein convergence order ~", order_mil)
print("State-dependent noise (GBM): Euler-Maruyama convergence order ~", order_em)

# ----------------------------------------------------------------------
# (3) Constant-noise check: dX = -theta X dt + sigma dW  (s'=0).
# Milstein must reduce EXACTLY to Euler-Maruyama on the same path.
# ----------------------------------------------------------------------
theta, sig_c = 1.0, 0.5
f_c = lambda x: -theta * x
s_c = lambda x: sig_c
sp_c = lambda x: 0.0
Nc = 500
dW_c = rng.normal(0.0, np.sqrt(T / Nc), size=Nc)
_, x_mil_c = milstein(f_c, s_c, sp_c, X0, T, Nc, dW=dW_c)
_, x_em_c = euler_maruyama(f_c, s_c, X0, T, Nc, dW=dW_c)
max_diff = np.max(np.abs(x_mil_c - x_em_c))
print("Constant noise: max |Milstein - EulerMaruyama| over trajectory =", max_diff)
print("Constant noise: Milstein reduces to Euler-Maruyama exactly? ", bool(max_diff == 0.0))

# Explanation of why the check confirms the result:
print("Explanation: Because the Milstein correction 0.5*s*s'*(dW^2-dt) is nonzero "
      "only when s'(X)!=0, it lifts the strong order from ~0.5 (Euler-Maruyama) to "
      "~1.0 for state-dependent noise yet leaves constant-noise trajectories "
      "bit-for-bit identical to Euler-Maruyama, confirming the term matters exactly "
      "when the noise depends on the state.")

# ----------------------------------------------------------------------
# Figure: gene trajectory + convergence log-log plot
# ----------------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.5))
ax1.plot(t_gene, traj, lw=0.8, color="C0")
ax1.set_xlabel("time"); ax1.set_ylabel("X (expression)")
ax1.set_title("Gene-circuit SDE via Milstein (6C.4)")

ax2.loglog(dts, err_em, "o-", label=f"Euler-Maruyama (order~{order_em:.2f})")
ax2.loglog(dts, err_mil, "s-", label=f"Milstein (order~{order_mil:.2f})")
ax2.loglog(dts, err_mil[0] * (dts / dts[0]), "k--", lw=0.8, label="slope 1 ref")
ax2.loglog(dts, err_em[0] * (dts / dts[0])**0.5, "k:", lw=0.8, label="slope 0.5 ref")
ax2.set_xlabel("dt"); ax2.set_ylabel("strong error at T")
ax2.set_title("Strong convergence (state-dependent noise)")
ax2.legend(fontsize=8)

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.3.1_s3.png", dpi=120)
