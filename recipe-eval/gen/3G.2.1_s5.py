import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Toggle-switch parameters (X and Y repress each other)
# ---------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12

# Right-hand sides of the ODEs
def dXdt(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def dYdt(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# ---------------------------------------------------------------
# Nullclines expressed as explicit functions
# ---------------------------------------------------------------
# Y-nullcline: dY/dt = 0  ->  Y as a function of X
def Y_nullcline(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# X-nullcline: dX/dt = 0  ->  X as a function of Y
def X_nullcline(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# ---------------------------------------------------------------
# Explicit cumulative trapezoidal integration (no library shortcut)
# U(s) = - integral of force ds, so steady states (force = 0) are extrema
# ---------------------------------------------------------------
def cumulative_potential(coord, force):
    U = np.zeros_like(coord)
    for i in range(1, len(coord)):
        # trapezoid on segment [i-1, i]; minus sign => potential
        U[i] = U[i - 1] - 0.5 * (force[i] + force[i - 1]) * (coord[i] - coord[i - 1])
    return U

# ---------------------------------------------------------------
# Path 1: walk ALONG the Y-nullcline, parameterized by X.
# Residual "force" driving the remaining variable is dX/dt.
# ---------------------------------------------------------------
X_path = np.linspace(0.0, (gX0 + gX1) / kX * 1.05, 4000)  # X from 0 to a bit above its max
Y_on_Ynull = Y_nullcline(X_path)                          # Y that keeps dY/dt = 0
force_X = dXdt(X_path, Y_on_Ynull)                        # net drift of X along this curve
U_Ynull = cumulative_potential(X_path, force_X)           # accumulated potential

# ---------------------------------------------------------------
# Path 2: walk ALONG the X-nullcline, parameterized by Y.
# Residual "force" driving the remaining variable is dY/dt.
# ---------------------------------------------------------------
Y_path = np.linspace(0.0, (gY0 + gY1) / kY * 1.05, 4000)  # Y from 0 to a bit above its max
X_on_Xnull = X_nullcline(Y_path)                          # X that keeps dX/dt = 0
force_Y = dYdt(X_on_Xnull, Y_path)                        # net drift of Y along this curve
U_Xnull = cumulative_potential(Y_path, force_Y)           # accumulated potential

# ---------------------------------------------------------------
# Locate steady states = zeros of the residual force = extrema of U.
# Detect sign changes along each path and refine linearly.
# ---------------------------------------------------------------
def find_extrema(coord, force, potential):
    roots = []
    for i in range(1, len(force)):
        if force[i - 1] == 0.0 or force[i - 1] * force[i] < 0.0:
            # linear interpolation of the zero crossing
            t = force[i - 1] / (force[i - 1] - force[i]) if force[i] != force[i - 1] else 0.0
            c = coord[i - 1] + t * (coord[i] - coord[i - 1])
            u = potential[i - 1] + t * (potential[i] - potential[i - 1])
            # slope of force tells stable (force decreasing -> min of U) vs saddle
            slope = (force[i] - force[i - 1]) / (coord[i] - coord[i - 1])
            kind = "stable (U minimum)" if slope < 0 else "saddle (U maximum)"
            roots.append((c, u, kind))
    return roots

extrema_Ynull = find_extrema(X_path, force_X, U_Ynull)
extrema_Xnull = find_extrema(Y_path, force_Y, U_Xnull)

# ---------------------------------------------------------------
# Report numerical results
# ---------------------------------------------------------------
print("=== Steady states located along the Y-nullcline (parameterized by X) ===")
for c, u, kind in extrema_Ynull:
    print(f"X = {c:10.4f}   Y = {Y_nullcline(c):10.4f}   U_Ynull = {u:12.4f}   -> {kind}")

print("\n=== Steady states located along the X-nullcline (parameterized by Y) ===")
for c, u, kind in extrema_Xnull:
    print(f"Y = {c:10.4f}   X = {X_nullcline(c):10.4f}   U_Xnull = {u:12.4f}   -> {kind}")

# Numerical comparison: the two paths give different potential values (flow is not a gradient)
print("\n=== Path-dependence check (potential values need not agree) ===")
for i, (c, u, kind) in enumerate(extrema_Ynull):
    print(f"Extremum {i}: Y-nullcline U = {u:12.4f}   [{kind}]")
for i, (c, u, kind) in enumerate(extrema_Xnull):
    print(f"Extremum {i}: X-nullcline U = {u:12.4f}   [{kind}]")

# ---------------------------------------------------------------
# Plots: accumulated potential along each nullcline, vs X and vs Y
# ---------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

# Y-nullcline potential vs X
ax[0, 0].plot(X_path, U_Ynull, 'b-')
for c, u, kind in extrema_Ynull:
    ax[0, 0].plot(c, u, 'ro' if "saddle" in kind else 'gs')
ax[0, 0].set_xlabel("X"); ax[0, 0].set_ylabel("Accumulated potential")
ax[0, 0].set_title("Along Y-nullcline  vs X")

# Y-nullcline potential vs Y (Y = Y_nullcline(X))
ax[0, 1].plot(Y_on_Ynull, U_Ynull, 'b-')
for c, u, kind in extrema_Ynull:
    ax[0, 1].plot(Y_nullcline(c), u, 'ro' if "saddle" in kind else 'gs')
ax[0, 1].set_xlabel("Y"); ax[0, 1].set_ylabel("Accumulated potential")
ax[0, 1].set_title("Along Y-nullcline  vs Y")

# X-nullcline potential vs X (X = X_nullcline(Y))
ax[1, 0].plot(X_on_Xnull, U_Xnull, 'm-')
for c, u, kind in extrema_Xnull:
    ax[1, 0].plot(X_nullcline(c), u, 'ro' if "saddle" in kind else 'gs')
ax[1, 0].set_xlabel("X"); ax[1, 0].set_ylabel("Accumulated potential")
ax[1, 0].set_title("Along X-nullcline  vs X")

# X-nullcline potential vs Y
ax[1, 1].plot(Y_path, U_Xnull, 'm-')
for c, u, kind in extrema_Xnull:
    ax[1, 1].plot(c, u, 'ro' if "saddle" in kind else 'gs')
ax[1, 1].set_xlabel("Y"); ax[1, 1].set_ylabel("Accumulated potential")
ax[1, 1].set_title("Along X-nullcline  vs Y")

plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3G.2.1_s5.png")

# One-sentence explanation of why the check confirms the result:
print("\nExplanation: Because the accumulated potential is defined as minus the integral "
      "of the residual drift, its stationary points occur exactly where that drift vanishes "
      "(the nullcline intersections), so minima at the stable states and a maximum at the "
      "saddle confirm we have correctly located the steady states.")
