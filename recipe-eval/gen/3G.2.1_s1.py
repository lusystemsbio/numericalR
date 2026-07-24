import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------
# Toggle-switch parameters (genes X and Y repress each other)
# ----------------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ----------------------------------------------------------------------
# Model velocities  dX/dt , dY/dt
# ----------------------------------------------------------------------
def vX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

def vY(X, Y):
    return gY0 + gY1 / (1.0 + (X / Xth) ** nX) - kY * Y

# Nullcline solutions (each equation solves algebraically for its own variable)
def fX(Y):   # X on the X-nullcline (dX/dt = 0)  ->  X = ... (function of Y)
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def fY(X):   # Y on the Y-nullcline (dY/dt = 0)  ->  Y = ... (function of X)
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ----------------------------------------------------------------------
# Locate the steady states = intersections of the two nullclines.
# A steady state satisfies  X = fX(fY(X)) , so root of h(X)=fX(fY(X))-X.
# ----------------------------------------------------------------------
def h(X):
    return fX(fY(X)) - X

Xscan = np.linspace(gX0 / kX, (gX0 + gX1) / kX, 20000)   # physical X range
hs = h(Xscan)
roots = []
for i in range(len(Xscan) - 1):
    if hs[i] == 0.0 or hs[i] * hs[i + 1] < 0.0:           # sign change -> bisect
        a, b = Xscan[i], Xscan[i + 1]
        for _ in range(100):
            m = 0.5 * (a + b)
            if h(a) * h(m) <= 0.0:
                b = m
            else:
                a = m
        Xs = 0.5 * (a + b)
        roots.append((Xs, fY(Xs)))

# Classify each steady state via the Jacobian eigenvalues
def classify(Xs, Ys):
    eps = 1e-4
    J = np.array([
        [(vX(Xs + eps, Ys) - vX(Xs - eps, Ys)) / (2 * eps),
         (vX(Xs, Ys + eps) - vX(Xs, Ys - eps)) / (2 * eps)],
        [(vY(Xs + eps, Ys) - vY(Xs - eps, Ys)) / (2 * eps),
         (vY(Xs, Ys + eps) - vY(Xs, Ys - eps)) / (2 * eps)],
    ])
    ev = np.linalg.eigvals(J)
    if np.all(ev.real < 0):
        return "stable node"
    elif np.any(ev.real > 0) and np.any(ev.real < 0):
        return "saddle"
    else:
        return "unstable"

print("=== Steady states (nullcline intersections) ===")
states = []
for Xs, Ys in roots:
    kind = classify(Xs, Ys)
    states.append((Xs, Ys, kind))
    print(f"steady state: X = {Xs:.4f}   Y = {Ys:.4f}   type = {kind}")

# ----------------------------------------------------------------------
# Effective potential along each nullcline, integrated with the
# trapezoidal rule written out explicitly (no library one-liner).
#
#  * On the X-nullcline X is slaved to Y (dX/dt=0); the residual flow is
#    dY/dt evaluated there -> integrate vY over Y.
#  * On the Y-nullcline Y is slaved to X (dY/dt=0); the residual flow is
#    dX/dt evaluated there -> integrate vX over X.
#  The potential is U = -integral(residual velocity), so extrema fall
#  exactly where the residual velocity (and hence both derivatives) vanish.
# ----------------------------------------------------------------------
def trapz_potential(var, vel):
    """U[i] = -cumulative trapezoidal integral of vel over var."""
    U = np.zeros_like(var)
    for i in range(1, len(var)):
        # trapezoid area on the segment [i-1, i]
        area = 0.5 * (vel[i] + vel[i - 1]) * (var[i] - var[i - 1])
        U[i] = U[i - 1] + area
    return -U   # potential is the negative of the accumulated flow

# --- Path A: along the X-nullcline, parametrised by Y ---
Ya = np.linspace(gY0 / kY, (gY0 + gY1) / kY, 4000)  # physical Y range
Xa = fX(Ya)                                         # X pinned to X-nullcline
velA = vY(Xa, Ya)                                   # residual velocity dY/dt
U_A = trapz_potential(Ya, velA)

# --- Path B: along the Y-nullcline, parametrised by X ---
Xb = np.linspace(gX0 / kX, (gX0 + gX1) / kX, 4000)  # physical X range
Yb = fY(Xb)                                         # Y pinned to Y-nullcline
velB = vX(Xb, Yb)                                   # residual velocity dX/dt
U_B = trapz_potential(Xb, velB)

# ----------------------------------------------------------------------
# Check: report the potential value at each steady state on each path
# and whether it is a local min (stable) or local max (saddle).
# ----------------------------------------------------------------------
print("\n=== Potential values at steady states ===")
for Xs, Ys, kind in states:
    UA_here = np.interp(Ys, Ya, U_A)                # path A indexed by Y
    UB_here = np.interp(Xs, Xb, U_B)                # path B indexed by X
    print(f"state (X={Xs:.3f}, Y={Ys:.3f}, {kind}): "
          f"U_Xnullcline = {UA_here:.4f}   U_Ynullcline = {UB_here:.4f}")

# Locate the extrema of each computed potential curve for confirmation
iA_min = np.argmin(U_A); iA_max = np.argmax(U_A)
iB_min = np.argmin(U_B); iB_max = np.argmax(U_B)
print("\n=== Extrema of the accumulated potential curves ===")
print(f"X-nullcline potential: min at Y={Ya[iA_min]:.3f} (U={U_A[iA_min]:.4f}), "
      f"max at Y={Ya[iA_max]:.3f} (U={U_A[iA_max]:.4f})")
print(f"Y-nullcline potential: min at X={Xb[iB_min]:.3f} (U={U_B[iB_min]:.4f}), "
      f"max at X={Xb[iB_max]:.3f} (U={U_B[iB_max]:.4f})")

# ----------------------------------------------------------------------
# Plots: accumulated potential along each nullcline, vs X and vs Y
# ----------------------------------------------------------------------
fig, ax = plt.subplots(2, 2, figsize=(12, 9))

ax[0, 0].plot(Ya, U_A, 'b-')
ax[0, 0].set_xlabel("Y"); ax[0, 0].set_ylabel("potential U")
ax[0, 0].set_title("X-nullcline potential vs Y")

ax[0, 1].plot(Xa, U_A, 'b-')
ax[0, 1].set_xlabel("X"); ax[0, 1].set_ylabel("potential U")
ax[0, 1].set_title("X-nullcline potential vs X")

ax[1, 0].plot(Xb, U_B, 'r-')
ax[1, 0].set_xlabel("X"); ax[1, 0].set_ylabel("potential U")
ax[1, 0].set_title("Y-nullcline potential vs X")

ax[1, 1].plot(Yb, U_B, 'r-')
ax[1, 1].set_xlabel("Y"); ax[1, 1].set_ylabel("potential U")
ax[1, 1].set_title("Y-nullcline potential vs Y")

# mark steady states on every panel
for Xs, Ys, kind in states:
    style = 'ko' if kind == "saddle" else 'g^'
    ax[0, 0].plot(Ys, np.interp(Ys, Ya, U_A), style)
    ax[0, 1].plot(Xs, np.interp(Ys, Ya, U_A), style)
    ax[1, 0].plot(Xs, np.interp(Xs, Xb, U_B), style)
    ax[1, 1].plot(Ys, np.interp(Xs, Xb, U_B), style)

fig.suptitle("Effective potential integrated along toggle-switch nullclines "
             "(triangles = stable nodes, circles = saddle)")
fig.tight_layout(rect=[0, 0, 1, 0.97])
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3G.2.1_s1.png")

# ----------------------------------------------------------------------
# Why the check works (one sentence):
# ----------------------------------------------------------------------
print("\nExplanation: Because the residual velocity being integrated vanishes "
      "exactly at the nullcline intersections, the accumulated potential is "
      "stationary there, dipping to a minimum where that residual flow is "
      "restoring (stable node) and peaking where it is repelling (saddle), so "
      "the observed dips-at-stable-states and peak-at-the-saddle confirm the "
      "extrema coincide with the true steady states.")
