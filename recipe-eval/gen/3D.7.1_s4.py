import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# --- Lynx-hare yearly data (Hudson's Bay Company, 1900-1920, thousands) ---
# Hare = prey N, Lynx = predator P.
data = np.array([
    [1900, 30.0,  4.0], [1901, 47.2,  6.1], [1902, 70.2,  9.8],
    [1903, 77.4, 35.2], [1904, 36.3, 59.4], [1905, 20.6, 41.7],
    [1906, 18.1, 19.0], [1907, 21.4, 13.0], [1908, 22.0,  8.3],
    [1909, 25.4,  9.1], [1910, 27.1,  7.4], [1911, 40.3,  8.0],
    [1912, 57.0, 12.3], [1913, 76.6, 19.5], [1914, 52.3, 45.7],
    [1915, 19.5, 51.1], [1916, 11.2, 29.7], [1917,  7.6, 15.8],
    [1918, 14.6,  9.7], [1919, 16.2, 10.1], [1920, 24.7,  8.6],
])
t = data[:, 0]      # years
N = data[:, 1]      # prey (hare)
P = data[:, 2]      # predator (lynx)

# --- Centered finite differences approximate the derivatives at interior pts ---
# dy/dt at i ~ (y[i+1]-y[i-1]) / (t[i+1]-t[i-1]); endpoints are dropped.
dNdt = (N[2:] - N[:-2]) / (t[2:] - t[:-2])
dPdt = (P[2:] - P[:-2]) / (t[2:] - t[:-2])
Ni, Pi = N[1:-1], P[1:-1]   # data values aligned with the derivatives

# --- The model is LINEAR in its parameters at each data point ---
# dN/dt = a*N - b*(N*P)   -> regress dNdt on columns [N, N*P] -> [a, -b]
# dP/dt = c*(N*P) - d*P   -> regress dPdt on columns [N*P, P] -> [c, -d]
X_prey = np.column_stack([Ni, Ni * Pi])      # design matrix, no intercept
X_pred = np.column_stack([Ni * Pi, Pi])

def ols(X, y):
    # Ordinary least squares via the normal equations: beta = (X^T X)^-1 X^T y
    XtX = X.T @ X
    Xty = X.T @ y
    return np.linalg.solve(XtX, Xty)         # closed-form, deterministic solve

beta_prey = ols(X_prey, dNdt)                # -> [a, -b]
beta_pred = ols(X_pred, dPdt)                # -> [c, -d]

a =  beta_prey[0]
b = -beta_prey[1]
c =  beta_pred[0]
d = -beta_pred[1]

print(f"a (prey growth rate)      = {a:.6f}")
print(f"b (predation rate)        = {b:.6f}")
print(f"c (predator growth rate)  = {c:.6f}")
print(f"d (predator death rate)   = {d:.6f}")

# --- Fitted trajectory: integrate the ODE from the first data point ---
def lv(y, tt, a, b, c, d):
    n, p = y
    return [n * (a - b * p), p * (c * n - d)]

t_dense = np.linspace(t[0], t[-1], 400)
sol = odeint(lv, [N[0], P[0]], t_dense, args=(a, b, c, d))
Nfit, Pfit = sol[:, 0], sol[:, 1]

# --- Check: the method returns parameters directly and deterministically ---
# Re-solving the same linear system must reproduce the parameters bit-for-bit,
# since OLS has a unique closed-form solution with no random initialization.
beta_prey_2 = ols(X_prey, dNdt)
beta_pred_2 = ols(X_pred, dPdt)
identical = (np.array_equal(beta_prey, beta_prey_2)
             and np.array_equal(beta_pred, beta_pred_2))
print(f"deterministic check: re-solve gives identical parameters = {identical}")
print(f"max param difference on re-solve = "
      f"{max(np.max(np.abs(beta_prey - beta_prey_2)), np.max(np.abs(beta_pred - beta_pred_2))):.3e}")
# One sentence: because OLS is a single closed-form linear solve with no random
# starting guess, re-running it on the same data must return the exact same
# numbers, which confirms the fit is produced directly and deterministically.

# --- Plot fitted trajectory over the data ---
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(t, N, "o", color="tab:green", label="Hare data (prey N)")
ax.plot(t, P, "s", color="tab:red", label="Lynx data (predator P)")
ax.plot(t_dense, Nfit, "-", color="tab:green", label="Fitted N(t)")
ax.plot(t_dense, Pfit, "-", color="tab:red", label="Fitted P(t)")
ax.set_xlabel("Year")
ax.set_ylabel("Population (thousands)")
ax.set_title("Lotka-Volterra fit by linear regression (finite-difference derivatives)")
ax.legend()
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.7.1_s4.png")
