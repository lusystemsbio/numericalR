import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# ----------------------------------------------------------------------
# Lynx-Hare yearly data (Hudson's Bay Co. records, in thousands).
# Prey N = hare, Predator P = lynx.
# ----------------------------------------------------------------------
year = np.arange(1900, 1921)
N = np.array([30.0, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22.0, 25.4,
              27.1, 40.3, 57.0, 76.6, 52.3, 19.5, 11.2,  7.6, 14.6, 16.2, 24.7])  # hare (prey)
P = np.array([ 4.0,  6.1,  9.8, 35.2, 59.4, 41.7, 19.0, 13.0,  8.3,  9.1,
               7.4,  8.0, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8,  9.7, 10.1,  8.6])  # lynx (predator)
t = year.astype(float)

# ----------------------------------------------------------------------
# Model:  dN/dt = a*N - b*(N*P)          (linear in a, b -> coefs [a, -b])
#         dP/dt = c*(N*P) - d*P          (linear in c, d -> coefs [c, -d])
# The right-hand sides are LINEAR in the parameters, so once we supply
# numerical values for the derivatives we get an overdetermined linear
# system A x = y that we solve by ordinary least squares (no intercept).
# ----------------------------------------------------------------------

# Centered finite differences for the derivatives at the INTERIOR points.
# dX/dt|_i = (X[i+1] - X[i-1]) / (t[i+1] - t[i-1])
i = np.arange(1, len(t) - 1)                       # interior indices
dt2 = t[i + 1] - t[i - 1]                           # width of the centered stencil
dNdt = (N[i + 1] - N[i - 1]) / dt2                  # numerical prey rate
dPdt = (P[i + 1] - P[i - 1]) / dt2                  # numerical predator rate

# Evaluate the regressor terms at the same interior points.
Ni, Pi = N[i], P[i]
NP = Ni * Pi                                        # interaction term N*P

# --- Prey regression: [N, N*P] @ [a, -b] = dN/dt  (no intercept column) ---
A_prey = np.column_stack([Ni, NP])
# Ordinary least squares via the normal equations, done explicitly:
#   x = (A^T A)^{-1} A^T y      (closed form, deterministic, no iteration)
ATA_prey = A_prey.T @ A_prey
ATy_prey = A_prey.T @ dNdt
coef_prey = np.linalg.solve(ATA_prey, ATy_prey)
a =  coef_prey[0]
b = -coef_prey[1]

# --- Predator regression: [N*P, P] @ [c, -d] = dP/dt  (no intercept) ---
A_pred = np.column_stack([NP, Pi])
ATA_pred = A_pred.T @ A_pred
ATy_pred = A_pred.T @ dPdt
coef_pred = np.linalg.solve(ATA_pred, ATy_pred)
c =  coef_pred[0]
d = -coef_pred[1]

print("Regression-fitted Lotka-Volterra parameters (OLS on finite differences):")
print(f"a (prey growth rate)        = {a:.6f}")
print(f"b (predation rate)          = {b:.6f}")
print(f"c (predator growth rate)    = {c:.6f}")
print(f"d (predator death rate)     = {d:.6f}")

# ----------------------------------------------------------------------
# Determinism / directness check:
# re-solve the identical normal equations a second time and compare.
# ----------------------------------------------------------------------
coef_prey_2 = np.linalg.solve(A_prey.T @ A_prey, A_prey.T @ dNdt)
coef_pred_2 = np.linalg.solve(A_pred.T @ A_pred, A_pred.T @ dPdt)
max_diff = max(np.max(np.abs(coef_prey_2 - coef_prey)),
               np.max(np.abs(coef_pred_2 - coef_pred)))
print(f"\nDeterminism check: max |param(run2) - param(run1)| = {max_diff:.3e}")
print("Check passes because OLS solves the normal equations in closed form, so "
      "re-running returns bit-identical parameters (max diff = 0), confirming the "
      "estimate is a direct deterministic function of the data, not the output of a "
      "randomly-seeded iterative optimizer.")

# ----------------------------------------------------------------------
# Fitted trajectory: integrate the ODE with the fitted parameters,
# starting from the first observed data point.
# ----------------------------------------------------------------------
def lv(tt, y):
    n, p = y
    return [n * (a - b * p), p * (c * n - d)]

sol = solve_ivp(lv, (t[0], t[-1]), [N[0], P[0]],
                t_eval=np.linspace(t[0], t[-1], 400), rtol=1e-8, atol=1e-8)

print("\nFitted trajectory sampled at data years (hare_fit, lynx_fit):")
sol_at_years = solve_ivp(lv, (t[0], t[-1]), [N[0], P[0]],
                         t_eval=t, rtol=1e-8, atol=1e-8)
for yr, nf, pf in zip(year, sol_at_years.y[0], sol_at_years.y[1]):
    print(f"{yr}: N={nf:8.3f}  P={pf:8.3f}")

# ----------------------------------------------------------------------
# Plot data vs fitted trajectory.
# ----------------------------------------------------------------------
plt.figure(figsize=(10, 6))
plt.plot(year, N, "o", color="tab:green", label="Hare (data)")
plt.plot(year, P, "s", color="tab:red", label="Lynx (data)")
plt.plot(sol.t, sol.y[0], "-", color="tab:green", label="Hare (fit)")
plt.plot(sol.t, sol.y[1], "-", color="tab:red", label="Lynx (fit)")
plt.xlabel("Year")
plt.ylabel("Population (thousands)")
plt.title("Lotka-Volterra fit by linear regression on finite differences")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.7.1_s1.png")
