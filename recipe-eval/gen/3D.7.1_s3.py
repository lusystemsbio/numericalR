import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# -----------------------------------------------------------------------------
# Lynx-hare yearly data (Hudson's Bay Company pelt counts, in thousands).
# Prey N = hare, predator P = lynx.
# -----------------------------------------------------------------------------
t = np.arange(1900, 1921, dtype=float)
N = np.array([30.0, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22.0, 25.4,
              27.1, 40.3, 57.0, 76.6, 52.3, 19.5, 11.2, 7.6, 14.6, 16.2, 24.7])
P = np.array([4.0, 6.1, 9.8, 35.2, 59.4, 41.7, 19.0, 13.0, 8.3, 9.1,
              7.4, 8.0, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7, 10.1, 8.6])

# -----------------------------------------------------------------------------
# The model is LINEAR IN ITS PARAMETERS:
#   dN/dt = a*N - b*(N*P)        (unknowns a, b enter linearly)
#   dP/dt = c*(N*P) - d*P        (unknowns c, d enter linearly)
# so if we replace the derivatives with numbers (finite differences) each data
# point becomes one linear equation in the parameters -> solve by OLS.
# -----------------------------------------------------------------------------

# Centered finite differences approximate dN/dt, dP/dt at the interior points i:
#   (f[i+1] - f[i-1]) / (t[i+1] - t[i-1])
i = np.arange(1, len(t) - 1)                      # interior indices
dt2 = t[i + 1] - t[i - 1]                          # spacing spanning 2 steps
dNdt = (N[i + 1] - N[i - 1]) / dt2                 # measured prey rate
dPdt = (P[i + 1] - P[i - 1]) / dt2                 # measured predator rate

# State values at the same interior points (right-hand-side "regressors").
Ni, Pi = N[i], P[i]

# --- Prey regression: dN/dt = [N, N*P] @ [a, -b]  (no intercept column) -------
X_prey = np.column_stack([Ni, Ni * Pi])            # design matrix, 2 columns
# Ordinary least squares via the normal equations: coef = (X^T X)^-1 X^T y
coef_prey = np.linalg.solve(X_prey.T @ X_prey, X_prey.T @ dNdt)
a = coef_prey[0]                                   # prey growth rate
b = -coef_prey[1]                                  # predation rate (sign flip)

# --- Predator regression: dP/dt = [N*P, P] @ [c, -d] --------------------------
X_pred = np.column_stack([Ni * Pi, Pi])            # design matrix, 2 columns
coef_pred = np.linalg.solve(X_pred.T @ X_pred, X_pred.T @ dPdt)
c = coef_pred[0]                                   # predator growth-on-prey rate
d = -coef_pred[1]                                  # predator death rate (sign flip)

print("Regression-fitted Lotka-Volterra parameters:")
print("a (prey growth)      =", a)
print("b (predation)        =", b)
print("c (predator growth)  =", c)
print("d (predator death)   =", d)

# -----------------------------------------------------------------------------
# Fitted trajectory: integrate the ODE forward from the first data point using
# the regression-recovered parameters, sampled at the data years.
# -----------------------------------------------------------------------------
def lv(tt, y):
    n, p = y
    return [n * (a - b * p), p * (c * n - d)]

sol = solve_ivp(lv, (t[0], t[-1]), [N[0], P[0]], t_eval=t,
                method="RK45", rtol=1e-8, atol=1e-8)
N_fit, P_fit = sol.y

print("\nFitted trajectory (year, hare_fit, lynx_fit):")
for yr, nf, pf in zip(t, N_fit, P_fit):
    print(f"{int(yr)}  {nf:8.3f}  {pf:8.3f}")

# -----------------------------------------------------------------------------
# CHECK: parameters are returned directly and deterministically.
# Re-solve the same OLS system and confirm bit-for-bit identical output.
# -----------------------------------------------------------------------------
coef_prey_again = np.linalg.solve(X_prey.T @ X_prey, X_prey.T @ dNdt)
coef_pred_again = np.linalg.solve(X_pred.T @ X_pred, X_pred.T @ dPdt)
identical = (np.array_equal(coef_prey, coef_prey_again)
             and np.array_equal(coef_pred, coef_pred_again))
print("\nDeterminism check (re-solve gives bit-identical params):", identical)
# One sentence: because OLS has a closed-form solution rather than an iterative
# search, re-running it reproduces the exact same numbers, which confirms the
# method returns the parameters directly and deterministically.

# -----------------------------------------------------------------------------
# Plot: data vs fitted trajectory.
# -----------------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(t, N, "o", color="tab:green", label="Hare data")
plt.plot(t, P, "s", color="tab:red", label="Lynx data")
plt.plot(t, N_fit, "-", color="tab:green", label="Hare fit")
plt.plot(t, P_fit, "-", color="tab:red", label="Lynx fit")
plt.xlabel("Year")
plt.ylabel("Population (thousands)")
plt.title("Lotka-Volterra fit by linear regression (lynx-hare data)")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.7.1_s3.png")
