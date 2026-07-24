import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Toggle-switch parameters
# ---------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# X-nullcline function fX(X,Y) = 0
def fX(X, Y):
    return gX0 + gX1 / (1.0 + (Y / Yth) ** nY) - kX * X

# Analytic partial derivatives of fX (the gradient components)
def dfX_dX(X, Y):
    return -kX
def dfX_dY(X, Y):
    # d/dY [ gX1 / (1 + (Y/Yth)^nY) ]
    r = (Y / Yth) ** nY
    dr_dY = nY * (Y / Yth) ** (nY - 1.0) * (1.0 / Yth)
    return -gX1 * dr_dY / (1.0 + r) ** 2

# ---------------------------------------------------------------
# Method 1: explicit arc-length continuation (predictor-corrector)
# State u = [X, Y]; follow the scalar constraint fX(u) = 0 by arc length.
# ---------------------------------------------------------------
def tangent(u, prev_t=None):
    # Gradient of fX
    g = np.array([dfX_dX(u[0], u[1]), dfX_dY(u[0], u[1])])
    # Tangent is perpendicular to the gradient (rotate by 90 deg)
    t = np.array([g[1], -g[0]])
    t = t / np.linalg.norm(t)
    # Keep the direction consistent with the previous step
    if prev_t is not None and np.dot(t, prev_t) < 0:
        t = -t
    return t

# Starting point on the curve: pick Y0 = 0, solve fX = 0 for X0 analytically
Y0 = 0.0
X0 = (gX0 + gX1 / (1.0 + (Y0 / Yth) ** nY)) / kX
u = np.array([X0, Y0])

ds = 4.0            # arc-length step
tol = 1e-12         # Newton tolerance
max_newton = 50     # max Newton iterations per corrector

# Force the initial tangent to move toward increasing Y
t = tangent(u, prev_t=np.array([0.0, 1.0]))

curveX, curveY = [u[0]], [u[1]]
for step in range(2000):
    # --- Predictor: step along the tangent by arc length ds ---
    u_pred = u + ds * t

    # --- Corrector: Newton on the augmented system ---
    #   F1 = fX(u)                     (stay on the nullcline)
    #   F2 = t . (u - u_pred)          (pseudo-arclength: stay near predictor)
    v = u_pred.copy()
    for _ in range(max_newton):
        F = np.array([fX(v[0], v[1]), np.dot(t, v - u_pred)])
        if np.linalg.norm(F) < tol:
            break
        J = np.array([[dfX_dX(v[0], v[1]), dfX_dY(v[0], v[1])],
                      [t[0], t[1]]])
        v = v - np.linalg.solve(J, F)

    # Update tangent using the new point, preserving direction
    t = tangent(v, prev_t=t)
    u = v
    curveX.append(u[0])
    curveY.append(u[1])

    # Stop when we have swept the region of interest
    if u[1] > 300.0 or u[0] < 0.0 or u[1] < -5.0:
        break

curveX = np.array(curveX)
curveY = np.array(curveY)
print(f"Continuation: number of traced points = {len(curveX)}")
print(f"Continuation: Y range traced = [{curveY.min():.4f}, {curveY.max():.4f}]")

# Residual of fX along the traced curve (should be ~0)
res_cont = np.max(np.abs(fX(curveX, curveY)))
print(f"Continuation: max |fX| along traced curve = {res_cont:.3e}")

# ---------------------------------------------------------------
# Method 2: direct separation (solve fX = 0 explicitly for X given Y)
# ---------------------------------------------------------------
def X_separation(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

Xsep = X_separation(curveY)
err_sep = np.max(np.abs(curveX - Xsep))
print(f"Separation vs continuation: max |X_cont - X_sep| = {err_sep:.3e}")

# ---------------------------------------------------------------
# Method 3: grid contour of fX at level 0
# ---------------------------------------------------------------
xg = np.linspace(0, 1200, 400)
yg = np.linspace(-5, 300, 400)
Xg, Yg = np.meshgrid(xg, yg)
Fg = fX(Xg, Yg)

cs = plt.contour(Xg, Yg, Fg, levels=[0.0])
segs = cs.allsegs[0]
contour_pts = np.vstack(segs)  # (N, 2) array of [X, Y] vertices
plt.close()

# Compare contour to the analytic curve at each contour vertex's Y
Xc, Yc = contour_pts[:, 0], contour_pts[:, 1]
err_contour = np.max(np.abs(Xc - X_separation(Yc)))
print(f"Contour vs separation: max |X_contour - X_sep| = {err_contour:.3e}")

# Cross-check contour against continuation by interpolating continuation X(Y)
order = np.argsort(curveY)
Xc_from_cont = np.interp(Yc, curveY[order], curveX[order])
err_contour_cont = np.max(np.abs(Xc - Xc_from_cont))
print(f"Contour vs continuation: max |X_contour - X_cont| = {err_contour_cont:.3e}")

# ---------------------------------------------------------------
# Phase-plane plot of the X-nullcline traced by continuation
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(curveX, curveY, 'b-', lw=2.5, label='X-nullcline (arc-length continuation)')
ax.plot(Xsep, curveY, 'r--', lw=1.2, label='separation X=(gX0+gX1/(1+(Y/Yth)^nY))/kX')
ax.plot(Xc, Yc, 'g.', ms=3, label='grid contour fX=0')
ax.plot(curveX[0], curveY[0], 'ko', ms=6, label='continuation start')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_title('Toggle switch: X-nullcline by arc-length continuation')
ax.legend(loc='upper right', fontsize=8)
ax.grid(True, alpha=0.3)
fig.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.6.1_s4.png")

# ---------------------------------------------------------------
# One-sentence explanation of why the agreement confirms the result.
# ---------------------------------------------------------------
print("Explanation: because all three methods must describe the same solution set "
      "{ (X,Y) : fX(X,Y)=0 }, their near-zero mutual discrepancies show the "
      "continuation curve lands on exactly the same nullcline, confirming it is correct.")
