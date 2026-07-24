import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Toggle-switch model parameters (genes X and Y mutually repress)
# dX/dt = gX0 + gX1/(1+(Y/Yth)^nY) - kX*X
# dY/dt = gY0 + gY1/(1+(X/Xth)^nX) - kY*Y
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Right-hand-side "flows" (production - degradation) for each variable
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def fY(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# Nullclines solved explicitly for the slaved variable:
# X-nullcline (dX/dt = 0): X as a function of Y
def X_on_Xnull(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# Y-nullcline (dY/dt = 0): Y as a function of X
def Y_on_Ynull(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ----------------------------------------------------------------------
# Explicit cumulative trapezoidal integration (no library one-liner)
# Returns the running integral of y w.r.t. t, starting at 0.
# ----------------------------------------------------------------------
def cumulative_trapz(y, t):
    U = np.zeros_like(y, dtype=float)
    for i in range(1, len(y)):
        # area of one trapezoid added to the accumulated value
        U[i] = U[i - 1] + 0.5 * (y[i] + y[i - 1]) * (t[i] - t[i - 1])
    return U

# ----------------------------------------------------------------------
# Effective potential ALONG THE X-NULLCLINE.
# On this curve dX/dt = 0, so the residual flow is the Y-flow.
# Parametrize by Y; the residual force is g(Y) = fY(X_nc(Y), Y).
# Potential U = -integral(residual force) dY  ->  extrema where force = 0.
# ----------------------------------------------------------------------
Y_par = np.linspace(0.0, 500.0, 4001)          # parameter along X-nullcline
X_par = X_on_Xnull(Y_par)                       # corresponding X on the curve
force_Xnull = fY(X_par, Y_par)                  # residual (Y) flow on the curve
U_Xnull = -cumulative_trapz(force_Xnull, Y_par) # accumulated effective potential

# ----------------------------------------------------------------------
# Effective potential ALONG THE Y-NULLCLINE.
# On this curve dY/dt = 0, so the residual flow is the X-flow.
# Parametrize by X; the residual force is h(X) = fX(X, Y_nc(X)).
# ----------------------------------------------------------------------
X_par2 = np.linspace(0.0, 700.0, 4001)          # parameter along Y-nullcline
Y_par2 = Y_on_Ynull(X_par2)                     # corresponding Y on the curve
force_Ynull = fX(X_par2, Y_par2)                # residual (X) flow on the curve
U_Ynull = -cumulative_trapz(force_Ynull, X_par2)# accumulated effective potential

# ----------------------------------------------------------------------
# Locate steady states as zeros of the residual flow (sign changes),
# refined by linear interpolation. These are the potential extrema.
# ----------------------------------------------------------------------
def find_zeros(force, param):
    zeros = []
    for i in range(1, len(force)):
        if force[i - 1] == 0.0:
            zeros.append(param[i - 1])
        elif force[i - 1] * force[i] < 0.0:   # sign change bracket
            # linear interpolation for the crossing location
            p = param[i - 1] - force[i - 1] * (param[i] - param[i - 1]) / (force[i] - force[i - 1])
            zeros.append(p)
    return zeros

# Steady states as (Y*) along the X-nullcline, and (X*) along the Y-nullcline
Y_states = find_zeros(force_Xnull, Y_par)
X_states = find_zeros(force_Ynull, X_par2)

# Full (X,Y) coordinates of each steady state (both nullclines agree there)
print("=== Steady states (extrema of the effective potential) ===")
ss_from_Xnull = [(X_on_Xnull(Ys), Ys) for Ys in Y_states]
ss_from_Ynull = [(Xs, Y_on_Ynull(Xs)) for Xs in X_states]
for k, (Xs, Ys) in enumerate(ss_from_Xnull):
    print(f"X-nullcline steady state {k}: X* = {Xs:.4f}, Y* = {Ys:.4f}")
for k, (Xs, Ys) in enumerate(ss_from_Ynull):
    print(f"Y-nullcline steady state {k}: X* = {Xs:.4f}, Y* = {Ys:.4f}")

# Classify each extremum as minimum (stable) or maximum (saddle) using the
# curvature of the accumulated potential at the nearest sampled index.
def classify(param_states, param_grid, U):
    labels = []
    for ps in param_states:
        i = int(np.argmin(np.abs(param_grid - ps)))
        i = max(1, min(len(U) - 2, i))
        curv = U[i - 1] - 2.0 * U[i] + U[i + 1]   # second-difference ~ curvature
        kind = "MINIMUM (stable)" if curv > 0 else "MAXIMUM (saddle)"
        labels.append((ps, U[i], kind))
    return labels

print("\n=== Potential along X-nullcline at its extrema ===")
for Ys, Uval, kind in classify(Y_states, Y_par, U_Xnull):
    print(f"Y* = {Ys:8.4f}  (X* = {X_on_Xnull(Ys):8.4f})  U = {Uval:12.4f}  -> {kind}")

print("\n=== Potential along Y-nullcline at its extrema ===")
for Xs, Uval, kind in classify(X_states, X_par2, U_Ynull):
    print(f"X* = {Xs:8.4f}  (Y* = {Y_on_Ynull(Xs):8.4f})  U = {Uval:12.4f}  -> {kind}")

# ----------------------------------------------------------------------
# Show the two paths do NOT give the same potential values at the states
# (the vector field is not a true gradient system).
# ----------------------------------------------------------------------
print("\n=== Path comparison (potential is path dependent) ===")
for Xs, Ys in ss_from_Xnull:
    iX = int(np.argmin(np.abs(Y_par - Ys)))
    iY = int(np.argmin(np.abs(X_par2 - Xs)))
    print(f"State near (X*={Xs:8.3f}, Y*={Ys:8.3f}): "
          f"U_Xnull = {U_Xnull[iX]:12.4f} ,  U_Ynull = {U_Ynull[iY]:12.4f}")

# ----------------------------------------------------------------------
# Plots: accumulated potential along each nullcline, versus X and versus Y
# ----------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# X-nullcline potential vs Y (its integration variable)
ax[0, 0].plot(Y_par, U_Xnull, color="C0")
ax[0, 0].scatter(Y_states, [U_Xnull[int(np.argmin(np.abs(Y_par - ys)))] for ys in Y_states],
                 color="k", zorder=5)
ax[0, 0].set_xlabel("Y"); ax[0, 0].set_ylabel("accumulated U")
ax[0, 0].set_title("Effective potential along X-nullcline  vs  Y")

# X-nullcline potential vs X (X also varies along the curve)
ax[0, 1].plot(X_par, U_Xnull, color="C0")
ax[0, 1].scatter([X_on_Xnull(ys) for ys in Y_states],
                 [U_Xnull[int(np.argmin(np.abs(Y_par - ys)))] for ys in Y_states],
                 color="k", zorder=5)
ax[0, 1].set_xlabel("X"); ax[0, 1].set_ylabel("accumulated U")
ax[0, 1].set_title("Effective potential along X-nullcline  vs  X")

# Y-nullcline potential vs X (its integration variable)
ax[1, 0].plot(X_par2, U_Ynull, color="C1")
ax[1, 0].scatter(X_states, [U_Ynull[int(np.argmin(np.abs(X_par2 - xs)))] for xs in X_states],
                 color="k", zorder=5)
ax[1, 0].set_xlabel("X"); ax[1, 0].set_ylabel("accumulated U")
ax[1, 0].set_title("Effective potential along Y-nullcline  vs  X")

# Y-nullcline potential vs Y (Y also varies along the curve)
ax[1, 1].plot(Y_par2, U_Ynull, color="C1")
ax[1, 1].scatter([Y_on_Ynull(xs) for xs in X_states],
                 [U_Ynull[int(np.argmin(np.abs(X_par2 - xs)))] for xs in X_states],
                 color="k", zorder=5)
ax[1, 1].set_xlabel("Y"); ax[1, 1].set_ylabel("accumulated U")
ax[1, 1].set_title("Effective potential along Y-nullcline  vs  Y")

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3G.2.1_s3.png")

# ----------------------------------------------------------------------
# One-sentence justification of the check:
# ----------------------------------------------------------------------
print("\nWhy the check confirms the result:")
print("Because the fixed points are exactly the zeros of the residual flow, "
      "they must sit at the extrema of the integrated potential, and the sign "
      "of the curvature (minima for the stable nodes, a maximum for the saddle) "
      "reproduces the known stability, so matching extrema to steady states "
      "verifies both their locations and their stability.")
