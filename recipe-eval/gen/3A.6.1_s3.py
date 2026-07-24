import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Toggle-switch parameters (mutual repression, Hill kinetics)
# ---------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12

# X-nullcline function: fX(X, Y) = 0 defines the curve we trace
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

# Partial derivatives (the Jacobian entries) of fX
def dfX_dX(X, Y):
    return -kX                                   # d/dX of fX
def dfX_dY(X, Y):
    # d/dY of gX1/(1+(Y/Yth)^nY)  via chain rule
    r = (Y / Yth) ** nY
    dr_dY = nY * (Y / Yth) ** (nY - 1) / Yth
    return -gX1 * dr_dY / (1.0 + r) ** 2

# ---------------------------------------------------------------
# Arc-length numerical continuation (predictor-corrector)
# state u = (X, Y); we follow fX(u) = 0 by stepping in arc length.
# ---------------------------------------------------------------
def continue_nullcline(X0, Y0, ds=2.0, n_steps=400, Ymax=350.0):
    pts = [(X0, Y0)]
    X, Y = X0, Y0
    prev_t = None
    for _ in range(n_steps):
        # --- tangent: perpendicular to gradient (fX_X, fX_Y) ---
        gx, gy = dfX_dX(X, Y), dfX_dY(X, Y)
        t = np.array([gy, -gx])          # tangent to the level set
        t /= np.hypot(*t)                # unit tangent
        # keep marching in a consistent direction (Y increasing here)
        if prev_t is not None and np.dot(t, prev_t) < 0:
            t = -t
        if prev_t is None and t[1] < 0:  # first step: force Y to increase
            t = -t
        prev_t = t

        # --- predictor: linear step of length ds along the tangent ---
        Xp, Yp = X + ds * t[0], Y + ds * t[1]

        # --- corrector: Newton on two equations ---
        #   g1 = fX(X,Y)                       (stay on the curve)
        #   g2 = (u - u_pred) . t = 0          (pseudo-arclength constraint)
        Xc, Yc = Xp, Yp
        for _ in range(50):
            g1 = fX(Xc, Yc)
            g2 = (Xc - Xp) * t[0] + (Yc - Yp) * t[1]
            J = np.array([[dfX_dX(Xc, Yc), dfX_dY(Xc, Yc)],
                          [t[0],           t[1]]])
            delta = np.linalg.solve(J, [-g1, -g2])
            Xc, Yc = Xc + delta[0], Yc + delta[1]
            if np.hypot(*delta) < 1e-12:
                break
        X, Y = Xc, Yc
        pts.append((X, Y))
        if Y > Ymax:
            break
    return np.array(pts)

# starting point: at Y0 = 0 the curve gives X exactly
Y0 = 0.0
X0 = (gX0 + gX1 / (1.0 + (Y0 / Yth) ** nY)) / kX
curve = continue_nullcline(X0, Y0)
Xc, Yc = curve[:, 0], curve[:, 1]

# ---------------------------------------------------------------
# Method 2 (separation): fX is linear in X, so solve X(Y) explicitly
# ---------------------------------------------------------------
def X_explicit(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

X_exact = X_explicit(Yc)

# ---------------------------------------------------------------
# Method 3 (grid contour): evaluate fX on a grid, take the 0-level
# ---------------------------------------------------------------
Xg = np.linspace(0, 1200, 400)
Yg = np.linspace(0, 350, 400)
XX, YY = np.meshgrid(Xg, Yg)
FF = fX(XX, YY)

# ---------------------------------------------------------------
# Check: max deviation of the traced curve from the explicit curve
# ---------------------------------------------------------------
max_dev_explicit = np.max(np.abs(Xc - X_exact))
# residual of fX along the traced curve (should be ~machine zero)
max_residual = np.max(np.abs(fX(Xc, Yc)))

print("Number of continuation points:", len(curve))
print("Start point (X0, Y0): {:.6f}, {:.6f}".format(X0, Y0))
print("End point   (X, Y):   {:.6f}, {:.6f}".format(Xc[-1], Yc[-1]))
print("Max |fX| residual along traced curve:", max_residual)
print("Max deviation traced vs explicit X(Y):", max_dev_explicit)

# ---------------------------------------------------------------
# Phase-plane plot
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
cs = ax.contour(XX, YY, FF, levels=[0], colors="0.6",
                linewidths=6, alpha=0.4)  # grid-contour nullcline
ax.plot(X_exact, Yc, "g--", lw=2, label="explicit (separation)")
ax.plot(Xc, Yc, "r.", ms=4, label="arc-length continuation")
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_title("Toggle switch: X-nullcline traced by continuation")
ax.legend()
ax.plot([], [], color="0.6", lw=6, alpha=0.4, label="grid contour")  # legend proxy
ax.legend()
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.6.1_s3.png")

# One-sentence explanation of why the check confirms the result:
print("Why the check confirms it: the continuation curve, the explicitly "
      "solved X(Y), and the grid 0-contour all describe the same set "
      "fX(X,Y)=0, so their agreement (deviation ~", 
      "{:.2e}".format(max_dev_explicit),
      ") means the predictor-corrector faithfully traced the true nullcline.")
