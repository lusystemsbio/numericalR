import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# ------------------------------------------------------------------
# Lynx-hare yearly data (classic Hudson Bay Company record, thousands)
# N = hare (prey), P = lynx (predator)
# ------------------------------------------------------------------
year = np.array([1900,1901,1902,1903,1904,1905,1906,1907,1908,1909,1910,
                 1911,1912,1913,1914,1915,1916,1917,1918,1919,1920], float)
N = np.array([30.0,47.2,70.2,77.4,36.3,20.6,18.1,21.4,22.0,25.4,27.1,
              40.3,57.0,76.6,52.3,19.5,11.2,7.6,14.6,16.2,24.7])   # hare (prey)
P = np.array([4.0,6.1,9.8,35.2,59.4,41.7,19.0,13.0,8.3,9.1,7.4,
              8.0,12.3,19.5,45.7,51.1,29.7,15.8,9.7,10.1,8.6])     # lynx (predator)
t = year - year[0]   # time in years starting at 0

# ------------------------------------------------------------------
# Step 1: replace time derivatives with CENTERED finite differences.
# For each interior point i:  dx/dt ~ (x[i+1]-x[i-1]) / (t[i+1]-t[i-1])
# ------------------------------------------------------------------
idx = np.arange(1, len(t) - 1)              # interior points only
dNdt = np.empty(len(idx))
dPdt = np.empty(len(idx))
for k, i in enumerate(idx):
    dt = t[i + 1] - t[i - 1]
    dNdt[k] = (N[i + 1] - N[i - 1]) / dt
    dPdt[k] = (P[i + 1] - P[i - 1]) / dt
Ni, Pi = N[idx], P[idx]                       # sample values at interior points

# ------------------------------------------------------------------
# Step 2: the model is LINEAR in the parameters, so build two linear systems.
#   dN/dt = a*N - b*(N*P)   -> regressors [N, N*P],  coeffs [a, -b]
#   dP/dt = c*(N*P) - d*P   -> regressors [N*P, P],  coeffs [c, -d]
# No intercept column is included (pure OLS through the origin).
# ------------------------------------------------------------------
X_prey = np.column_stack([Ni, Ni * Pi])       # design matrix for prey equation
y_prey = dNdt
X_pred = np.column_stack([Ni * Pi, Pi])       # design matrix for predator equation
y_pred = dPdt

# ------------------------------------------------------------------
# Step 3: solve each OLS problem explicitly via the normal equations
#   beta = (X^T X)^{-1} X^T y   (closed form, no iteration, no intercept)
# ------------------------------------------------------------------
def ols_normal_equations(X, y):
    XtX = X.T @ X                             # Gram matrix
    Xty = X.T @ y                             # projection of y onto columns
    return np.linalg.solve(XtX, Xty)          # direct linear solve

beta_prey = ols_normal_equations(X_prey, y_prey)   # [a, -b]
beta_pred = ols_normal_equations(X_pred, y_pred)   # [c, -d]

a = beta_prey[0]
b = -beta_prey[1]
c = beta_pred[0]
d = -beta_pred[1]

print("Regression-fitted Lotka-Volterra parameters:")
print(f"a (prey growth rate)       = {a:.6f}")
print(f"b (predation rate)         = {b:.6f}")
print(f"c (predator growth rate)   = {c:.6f}")
print(f"d (predator death rate)    = {d:.6f}")

# ------------------------------------------------------------------
# Step 4: fitted trajectory -- integrate the ODE with the fitted params,
# starting from the first data point.
# ------------------------------------------------------------------
def lv(tt, y):
    n, p = y
    return [n * (a - b * p), p * (c * n - d)]

sol = solve_ivp(lv, (t[0], t[-1]), [N[0], P[0]],
                t_eval=np.linspace(t[0], t[-1], 400), rtol=1e-8, atol=1e-8)

print("\nFitted trajectory sampled at data years:")
fit_at_years = solve_ivp(lv, (t[0], t[-1]), [N[0], P[0]],
                         t_eval=t, rtol=1e-8, atol=1e-8)
for yr, nh, pl in zip(year, fit_at_years.y[0], fit_at_years.y[1]):
    print(f"year {int(yr)}: hare_fit = {nh:8.3f}, lynx_fit = {pl:8.3f}")

# ------------------------------------------------------------------
# Step 5: determinism check -- solve the SAME OLS a second time and
# confirm the parameters are returned directly and bit-for-bit identical.
# ------------------------------------------------------------------
beta_prey_2 = ols_normal_equations(X_prey, y_prey)
beta_pred_2 = ols_normal_equations(X_pred, y_pred)
max_diff = max(np.max(np.abs(beta_prey - beta_prey_2)),
               np.max(np.abs(beta_pred - beta_pred_2)))
print("\nDeterminism check:")
print(f"max parameter difference between two solves = {max_diff:.3e}")
print(f"deterministic (difference exactly zero): {max_diff == 0.0}")
# One sentence: because OLS has a closed-form solution with no random
# initialization or iteration, re-solving reproduces the parameters
# exactly (zero difference), confirming the fit is direct and deterministic
# rather than the product of a stochastic/iterative optimizer.

# ------------------------------------------------------------------
# Plot data + fitted trajectory
# ------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(year, N, "o", color="tab:green", label="hare data")
ax.plot(year, P, "s", color="tab:red", label="lynx data")
ax.plot(sol.t + year[0], sol.y[0], "-", color="tab:green", label="hare fit")
ax.plot(sol.t + year[0], sol.y[1], "-", color="tab:red", label="lynx fit")
ax.set_xlabel("year")
ax.set_ylabel("population (thousands)")
ax.set_title("Lotka-Volterra fit by finite-difference linear regression")
ax.legend()
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.7.1_s5.png")
