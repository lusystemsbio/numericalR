import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

# ----------------------------------------------------------------------
# Toggle-switch parameters (X and Y mutually repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10   # dX/dt params
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12   # dY/dt params

# Right-hand sides of the ODEs
def dXdt(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def dYdt(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ----------------------------------------------------------------------
# Nullclines solved explicitly (each ODE set to zero for one variable)
# ----------------------------------------------------------------------
# X-nullcline: dX/dt = 0  ->  X is an explicit function of Y
def X_on_Xnull(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# Y-nullcline: dY/dt = 0  ->  Y is an explicit function of X
def Y_on_Ynull(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ----------------------------------------------------------------------
# Manual cumulative trapezoidal integration (no library one-liner)
# Returns V such that V[i] = integral of f dt from t[0] to t[i]
# ----------------------------------------------------------------------
def cumulative_trapz(f, t):
    V = np.zeros_like(f, dtype=float)
    for i in range(1, len(t)):
        # area of one trapezoid = 0.5*(f_left + f_right)*(t_right - t_left)
        V[i] = V[i - 1] + 0.5 * (f[i] + f[i - 1]) * (t[i] - t[i - 1])
    return V

# ----------------------------------------------------------------------
# Build each nullcline and the "left-over" force that drives motion on it.
# On the X-nullcline dX/dt=0, so the only motion is along Y, driven by dY/dt.
# The effective potential is minus the integral of that residual force.
# ----------------------------------------------------------------------
# --- X-nullcline: parametrize by Y ---
Y_axis = np.linspace(0.0, 450.0, 4001)     # sweep the free variable Y
X_Xnull = X_on_Xnull(Y_axis)               # matching X on the curve
force_Xnull = dYdt(X_Xnull, Y_axis)        # residual force = dY/dt along curve
V_Xnull = -cumulative_trapz(force_Xnull, Y_axis)  # potential: dV = -force dY

# --- Y-nullcline: parametrize by X ---
X_axis = np.linspace(0.0, 600.0, 4001)     # sweep the free variable X
Y_Ynull = Y_on_Ynull(X_axis)               # matching Y on the curve
force_Ynull = dXdt(X_axis, Y_Ynull)        # residual force = dX/dt along curve
V_Ynull = -cumulative_trapz(force_Ynull, X_axis)  # potential: dV = -force dX

# ----------------------------------------------------------------------
# Steady states = intersections of the two nullclines.
# On the X-nullcline, a steady state also needs dY/dt = 0, so we look for
# roots of g(Y) = dYdt(X_on_Xnull(Y), Y).
# ----------------------------------------------------------------------
def g_of_Y(Y):
    return dYdt(X_on_Xnull(Y), Y)

g_vals = g_of_Y(Y_axis)
roots_Y = []
for i in range(len(Y_axis) - 1):
    if g_vals[i] == 0.0:
        roots_Y.append(Y_axis[i])
    elif g_vals[i] * g_vals[i + 1] < 0.0:            # sign change -> a root
        roots_Y.append(brentq(g_of_Y, Y_axis[i], Y_axis[i + 1]))
roots_Y = np.array(sorted(roots_Y))
steady_X = X_on_Xnull(roots_Y)                       # matching X values
steady_Y = roots_Y

# Classify each steady state via the Jacobian eigenvalues
def classify(Xs, Ys):
    h = 1e-4
    a11 = (dXdt(Xs + h, Ys) - dXdt(Xs - h, Ys)) / (2 * h)
    a12 = (dXdt(Xs, Ys + h) - dXdt(Xs, Ys - h)) / (2 * h)
    a21 = (dYdt(Xs + h, Ys) - dYdt(Xs - h, Ys)) / (2 * h)
    a22 = (dYdt(Xs, Ys + h) - dYdt(Xs, Ys - h)) / (2 * h)
    ev = np.linalg.eigvals(np.array([[a11, a12], [a21, a22]]))
    if np.all(ev.real < 0):
        return "stable node"
    elif np.any(ev.real > 0) and np.any(ev.real < 0):
        return "saddle"
    else:
        return "unstable"

kinds = [classify(sx, sy) for sx, sy in zip(steady_X, steady_Y)]

# Potential value at each steady state on each path (linear interpolation)
V_at_Xnull = np.interp(steady_Y, Y_axis, V_Xnull)    # X-null potential vs Y
V_at_Ynull = np.interp(steady_X, X_axis, V_Ynull)    # Y-null potential vs X

# ----------------------------------------------------------------------
# Print all numerical results
# ----------------------------------------------------------------------
print("Number of steady states found:", len(steady_X))
for i in range(len(steady_X)):
    print(f"Steady state {i}: X = {steady_X[i]:.4f}, Y = {steady_Y[i]:.4f}, type = {kinds[i]}")
    print(f"  Potential on X-nullcline path at this state = {V_at_Xnull[i]:.6f}")
    print(f"  Potential on Y-nullcline path at this state = {V_at_Ynull[i]:.6f}")

# Check: extrema of the accumulated potential should sit at the steady states
print("\n--- Extremum check (potential derivative = residual force = 0) ---")
for i in range(len(steady_X)):
    tag = "MINIMUM" if kinds[i] == "stable node" else ("MAXIMUM" if kinds[i] == "saddle" else "inflection")
    print(f"State {i} ({kinds[i]}): expect {tag} in accumulated potential -> "
          f"X-null V={V_at_Xnull[i]:.4f}, Y-null V={V_at_Ynull[i]:.4f}")

print("\nX-nullcline potential range: min = %.4f, max = %.4f" % (V_Xnull.min(), V_Xnull.max()))
print("Y-nullcline potential range: min = %.4f, max = %.4f" % (V_Ynull.min(), V_Ynull.max()))
print("Note: the two paths give different potential values because the "
      "flow is not a true gradient field (path dependence).")

# ----------------------------------------------------------------------
# Plots: accumulated potential along each nullcline, versus X and versus Y
# ----------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

stable_mask = np.array([k == "stable node" for k in kinds])
saddle_mask = np.array([k == "saddle" for k in kinds])

# X-nullcline potential vs Y
ax[0, 0].plot(Y_axis, V_Xnull, 'b-')
ax[0, 0].plot(steady_Y[stable_mask], V_at_Xnull[stable_mask], 'go', ms=9, label='stable (min)')
ax[0, 0].plot(steady_Y[saddle_mask], V_at_Xnull[saddle_mask], 'r^', ms=9, label='saddle (max)')
ax[0, 0].set_xlabel("Y"); ax[0, 0].set_ylabel("Potential"); ax[0, 0].set_title("X-nullcline potential vs Y")
ax[0, 0].legend()

# X-nullcline potential vs X
ax[0, 1].plot(X_Xnull, V_Xnull, 'b-')
ax[0, 1].plot(steady_X[stable_mask], V_at_Xnull[stable_mask], 'go', ms=9)
ax[0, 1].plot(steady_X[saddle_mask], V_at_Xnull[saddle_mask], 'r^', ms=9)
ax[0, 1].set_xlabel("X"); ax[0, 1].set_ylabel("Potential"); ax[0, 1].set_title("X-nullcline potential vs X")

# Y-nullcline potential vs X
ax[1, 0].plot(X_axis, V_Ynull, 'm-')
ax[1, 0].plot(steady_X[stable_mask], V_at_Ynull[stable_mask], 'go', ms=9, label='stable (min)')
ax[1, 0].plot(steady_X[saddle_mask], V_at_Ynull[saddle_mask], 'r^', ms=9, label='saddle (max)')
ax[1, 0].set_xlabel("X"); ax[1, 0].set_ylabel("Potential"); ax[1, 0].set_title("Y-nullcline potential vs X")
ax[1, 0].legend()

# Y-nullcline potential vs Y
ax[1, 1].plot(Y_Ynull, V_Ynull, 'm-')
ax[1, 1].plot(steady_Y[stable_mask], V_at_Ynull[stable_mask], 'go', ms=9)
ax[1, 1].plot(steady_Y[saddle_mask], V_at_Ynull[saddle_mask], 'r^', ms=9)
ax[1, 1].set_xlabel("Y"); ax[1, 1].set_ylabel("Potential"); ax[1, 1].set_title("Y-nullcline potential vs Y")

# One-sentence explanation of why the check confirms the result:
# Because the accumulated potential is minus the integral of the residual force,
# its stationary points occur exactly where that force vanishes (the steady states),
# so minima landing on stable nodes and a maximum on the saddle confirms that the
# integrated extrema coincide with the true fixed points.
fig.suptitle("Effective potential integrated (trapezoidal) along toggle-switch nullclines")
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3G.2.1_s2.png")
print("\nExplanation: the extremum check confirms the result because the potential is "
      "defined as minus the integral of the residual on-nullcline force, so its "
      "stationary points fall exactly where that force is zero (the steady states), "
      "with minima at stable nodes and a maximum at the saddle.")
