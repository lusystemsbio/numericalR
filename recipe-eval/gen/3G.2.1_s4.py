import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Toggle-switch parameters (X and Y mutually repress each other)
# ---------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# Right-hand sides of the ODEs
def dXdt(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def dYdt(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ---------------------------------------------------------------
# Nullcline expressions (solve each dInterval=0 for its own variable)
# ---------------------------------------------------------------
# Y-nullcline: dY/dt = 0  ->  Y as an explicit function of X
def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# X-nullcline: dX/dt = 0  ->  X as an explicit function of Y
def X_of_Y(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# ---------------------------------------------------------------
# Effective potential by explicit trapezoidal accumulation.
# Along a nullcline the "own" component is exactly zero, so the
# residual dynamics is the OTHER component; the potential is minus
# the running integral of that residual force along the path.
# ---------------------------------------------------------------
def trapz_potential(t, force):
    """U(t) = -integral_0^t force dt', accumulated by the trapezoid rule."""
    U = np.zeros_like(t)
    for i in range(1, len(t)):
        # trapezoid area of the force over [t[i-1], t[i]]
        area = 0.5 * (force[i] + force[i - 1]) * (t[i] - t[i - 1])
        U[i] = U[i - 1] - area          # potential is minus the work done by the force
    return U

# ---- Path 1: walk along the Y-nullcline, parametrised by X -----
X_path = np.linspace(0.0, 700.0, 4000)      # X sweep
Y_on_Ynull = Y_of_X(X_path)                 # Y sits on its own nullcline
res_X = dXdt(X_path, Y_on_Ynull)            # residual (non-zero) flow is dX/dt
U_alongX = trapz_potential(X_path, res_X)   # potential accumulated vs X

# ---- Path 2: walk along the X-nullcline, parametrised by Y -----
Y_path = np.linspace(0.0, 600.0, 4000)      # Y sweep
X_on_Xnull = X_of_Y(Y_path)                 # X sits on its own nullcline
res_Y = dYdt(X_on_Xnull, Y_path)            # residual (non-zero) flow is dY/dt
U_alongY = trapz_potential(Y_path, res_Y)   # potential accumulated vs Y

# ---------------------------------------------------------------
# Locate steady states as the extrema of the potential:
# they are exactly the zeros of the residual force (sign changes).
# ---------------------------------------------------------------
def find_zeros(t, f):
    """Linearly-interpolated sign-change crossings of f(t)."""
    zeros = []
    for i in range(1, len(f)):
        if f[i - 1] == 0.0:
            zeros.append(t[i - 1])
        elif f[i - 1] * f[i] < 0.0:
            # linear interpolation to the crossing point
            zeros.append(t[i - 1] - f[i - 1] * (t[i] - t[i - 1]) / (f[i] - f[i - 1]))
    return np.array(zeros)

X_states = find_zeros(X_path, res_X)   # steady-state X values (extrema of U vs X)
Y_states = find_zeros(Y_path, res_Y)   # steady-state Y values (extrema of U vs Y)

# Classify each extremum: a stable node is a minimum (residual force goes +->-),
# the saddle is a maximum (residual force goes -->+).
def classify(t, f, roots):
    labels = []
    for r in roots:
        j = np.searchsorted(t, r)
        slope = (f[min(j, len(f) - 1)] - f[max(j - 1, 0)])
        labels.append("stable (min)" if slope < 0 else "saddle (max)")
    return labels

X_labels = classify(X_path, res_X, X_states)
Y_labels = classify(Y_path, res_Y, Y_states)

# potential value at each located extremum (interpolated)
UX_at = np.interp(X_states, X_path, U_alongX)
UY_at = np.interp(Y_states, Y_path, U_alongY)

# ---------------------------------------------------------------
# Report every numerical result
# ---------------------------------------------------------------
print("=== Steady states along the Y-nullcline (potential vs X) ===")
for x, lab, u in zip(X_states, X_labels, UX_at):
    y = Y_of_X(x)
    print(f"X = {x:10.4f}   Y = {y:10.4f}   type = {lab:14s}   U = {u:12.4f}")

print("\n=== Steady states along the X-nullcline (potential vs Y) ===")
for y, lab, u in zip(Y_states, Y_labels, UY_at):
    x = X_of_Y(y)
    print(f"Y = {y:10.4f}   X = {x:10.4f}   type = {lab:14s}   U = {u:12.4f}")

print("\n=== Potential-value comparison at the extrema (paths need NOT agree) ===")
print(f"U(vs X) at extrema: {np.round(UX_at, 4)}")
print(f"U(vs Y) at extrema: {np.round(UY_at, 4)}")
print(f"Number of steady states found (X path): {len(X_states)}")
print(f"Number of steady states found (Y path): {len(Y_states)}")

# ---------------------------------------------------------------
# Plots of the accumulated potential along each nullcline
# ---------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))

ax[0].plot(X_path, U_alongX, color="tab:blue", lw=1.8)
for x, lab in zip(X_states, X_labels):
    u = np.interp(x, X_path, U_alongX)
    c = "green" if "stable" in lab else "red"
    ax[0].plot(x, u, "o", color=c, ms=8)
    ax[0].annotate(lab, (x, u), textcoords="offset points", xytext=(6, 8), fontsize=8)
ax[0].set_xlabel("X")
ax[0].set_ylabel("accumulated potential  U(X)")
ax[0].set_title("Potential along Y-nullcline (residual = dX/dt), vs X")
ax[0].grid(alpha=0.3)

ax[1].plot(Y_path, U_alongY, color="tab:purple", lw=1.8)
for y, lab in zip(Y_states, Y_labels):
    u = np.interp(y, Y_path, U_alongY)
    c = "green" if "stable" in lab else "red"
    ax[1].plot(y, u, "o", color=c, ms=8)
    ax[1].annotate(lab, (y, u), textcoords="offset points", xytext=(6, 8), fontsize=8)
ax[1].set_xlabel("Y")
ax[1].set_ylabel("accumulated potential  U(Y)")
ax[1].set_title("Potential along X-nullcline (residual = dY/dt), vs Y")
ax[1].grid(alpha=0.3)

fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3G.2.1_s4.png", dpi=130)

# One-sentence explanation of why the check confirms the result:
print("\nWhy the check confirms the result: because a steady state is exactly a zero "
      "of the residual force, the trapezoidal potential must be stationary there, "
      "dipping to a minimum where the force restores (stable node) and peaking where "
      "it repels (saddle); the two paths giving different potential VALUES confirms "
      "the flow is not a true gradient, so the potential is only an effective, "
      "path-dependent construct that still correctly pinpoints the extrema.")
