import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Lynx-hare yearly data (Hudson's Bay Company, 1900-1920), thousands of pelts.
# N = hare (prey), P = lynx (predator)
# ----------------------------------------------------------------------
year = np.arange(1900, 1921)
N = np.array([30.0, 47.2, 70.2, 77.4, 36.3, 20.6, 18.1, 21.4, 22.0, 25.4,
              27.1, 40.3, 57.0, 76.6, 52.3, 19.5, 11.2, 7.6, 14.6, 16.2, 24.7])   # prey (hare)
P = np.array([4.0, 6.1, 9.8, 35.2, 59.4, 41.7, 19.0, 13.0, 8.3, 9.1,
              7.4, 8.0, 12.3, 19.5, 45.7, 51.1, 29.7, 15.8, 9.7, 10.1, 8.6])       # predator (lynx)
t = year.astype(float)

# ----------------------------------------------------------------------
# The Lotka-Volterra model is LINEAR IN ITS PARAMETERS:
#   dN/dt = a*N - b*(N*P)      -> unknowns (a, b), regressors (N, N*P)
#   dP/dt = c*(N*P) - d*P      -> unknowns (c, d), regressors (N*P, P)
# So each rate equation becomes an ordinary-least-squares problem with no intercept.
# ----------------------------------------------------------------------

# Step 1: replace the time derivatives with CENTERED finite differences at
# each interior data point:  df/dt|_i ~ (f_{i+1} - f_{i-1}) / (t_{i+1} - t_{i-1})
dt2 = t[2:] - t[:-2]                     # spacing across the centered stencil
dNdt = (N[2:] - N[:-2]) / dt2            # estimated dN/dt at points 1..n-2
dPdt = (P[2:] - P[:-2]) / dt2            # estimated dP/dt at points 1..n-2

# Interior samples of the state (aligned with the centered differences)
Ni = N[1:-1]
Pi = P[1:-1]
NPi = Ni * Pi                            # the bilinear interaction term

# Step 2: build the two design matrices (columns are the regressors, no intercept)
X_prey = np.column_stack([Ni, NPi])     # [N, N*P]  -> coeffs are (a, -b)
X_pred = np.column_stack([NPi, Pi])     # [N*P, P]  -> coeffs are (c, -d)

# Step 3: solve each linear system by ordinary least squares (rcond=None -> exact OLS)
coef_prey, *_ = np.linalg.lstsq(X_prey, dNdt, rcond=None)
coef_pred, *_ = np.linalg.lstsq(X_pred, dPdt, rcond=None)

# Step 4: read the physical parameters back out of the fitted coefficients
a = coef_prey[0]
b = -coef_prey[1]
c = coef_pred[0]
d = -coef_pred[1]

print(f"Fitted a (prey growth rate)     = {a:.6f}")
print(f"Fitted b (predation rate)       = {b:.6f}")
print(f"Fitted c (predator growth rate) = {c:.6f}")
print(f"Fitted d (predator death rate)  = {d:.6f}")

# ----------------------------------------------------------------------
# Determinism / directness check: re-solve via the closed-form normal
# equations  beta = (X^T X)^{-1} X^T y  and confirm it reproduces lstsq
# exactly (to machine precision, no iteration, no random seed).
# ----------------------------------------------------------------------
beta_prey_normal = np.linalg.solve(X_prey.T @ X_prey, X_prey.T @ dNdt)
beta_pred_normal = np.linalg.solve(X_pred.T @ X_pred, X_pred.T @ dPdt)
max_diff = max(np.max(np.abs(beta_prey_normal - coef_prey)),
               np.max(np.abs(beta_pred_normal - coef_pred)))
print(f"Max |lstsq - normal-equations| coefficient difference = {max_diff:.3e}")
# This confirms the result because the normal equations give THE unique OLS
# solution in closed form, so agreement to machine precision proves the
# parameters were returned directly and deterministically, not iteratively.

# ----------------------------------------------------------------------
# Fitted trajectory: integrate the ODE with the regressed parameters from the
# first data point, using a simple fixed-step RK4 (no external solver needed).
# ----------------------------------------------------------------------
def lv(y):
    n, p = y
    return np.array([n * (a - b * p), p * (c * n - d)])

def integrate(y0, t_start, t_end, h=0.01):
    ts = [t_start]
    ys = [np.array(y0, dtype=float)]
    tc, yc = t_start, np.array(y0, dtype=float)
    while tc < t_end - 1e-9:
        step = min(h, t_end - tc)
        k1 = lv(yc)
        k2 = lv(yc + 0.5 * step * k1)
        k3 = lv(yc + 0.5 * step * k2)
        k4 = lv(yc + step * k3)
        yc = yc + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        tc += step
        ts.append(tc)
        ys.append(yc)
    return np.array(ts), np.array(ys)

t_fit, y_fit = integrate([N[0], P[0]], t[0], t[-1], h=0.01)

# Report the fitted trajectory sampled at the data years
print("\nFitted trajectory at each data year (year, hare_fit, lynx_fit):")
for yr in year:
    idx = np.argmin(np.abs(t_fit - yr))
    print(f"{yr}  N_fit={y_fit[idx,0]:.3f}  P_fit={y_fit[idx,1]:.3f}")

# ----------------------------------------------------------------------
# Plot data vs fitted trajectory
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(year, N, "o", color="tab:green", label="Hare data (prey N)")
ax.plot(year, P, "s", color="tab:red", label="Lynx data (predator P)")
ax.plot(t_fit, y_fit[:, 0], "-", color="tab:green", label="Hare fit")
ax.plot(t_fit, y_fit[:, 1], "-", color="tab:red", label="Lynx fit")
ax.set_xlabel("Year")
ax.set_ylabel("Population (thousands)")
ax.set_title("Lotka-Volterra fit by linear regression (finite-difference derivatives, OLS)")
ax.legend()
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3D.7.1_s2.png")
print("\nSaved figure to 3D.7.1_s2.png")
