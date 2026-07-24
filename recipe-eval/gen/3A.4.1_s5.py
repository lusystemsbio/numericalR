import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Toggle switch parameters ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.10
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ---- Closed-form nullclines via separation of variables ----
# fX(X,Y) = gX0 + gX1/(1+(Y/Yth)^nY) - kX*X = 0
#   -> solve for X given Y:  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX
def X_nullcline(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# fY(X,Y) = gY0 + gY1/(1+(X/Xth)^nX) - kY*Y = 0
#   -> solve for Y given X:  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY
def Y_nullcline(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ---- Sweep the free variable for each nullcline ----
# X-nullcline: sweep Y, get corresponding X  -> curve of points (X(Y), Y)
Y_sweep = np.linspace(0.0, 600.0, 2000)
X_of_Y  = X_nullcline(Y_sweep)

# Y-nullcline: sweep X, get corresponding Y  -> curve of points (X, Y(X))
X_sweep = np.linspace(0.0, 600.0, 2000)
Y_of_X  = Y_nullcline(X_sweep)

# ---- Find intersections (steady states) ----
# At a crossing both hold: X = X_nullcline(Y) and Y = Y_nullcline(X).
# Substitute -> define residual g(X) = X_nullcline(Y_nullcline(X)) - X = 0.
# Locate roots by scanning for sign changes, then refine by bisection.
def g(X):
    return X_nullcline(Y_nullcline(X)) - X

Xscan = np.linspace(1.0, 600.0, 4000)
gv = g(Xscan)
roots = []
for i in range(len(Xscan) - 1):
    if gv[i] == 0.0:
        roots.append(Xscan[i])
    elif gv[i] * gv[i + 1] < 0.0:
        a, b = Xscan[i], Xscan[i + 1]
        for _ in range(100):  # bisection refinement
            m = 0.5 * (a + b)
            if g(a) * g(m) <= 0.0:
                b = m
            else:
                a = m
        roots.append(0.5 * (a + b))

# De-duplicate near-identical roots
roots_clean = []
for r in roots:
    if not any(abs(r - s) < 1e-3 for s in roots_clean):
        roots_clean.append(r)

steady_states = [(r, Y_nullcline(r)) for r in roots_clean]

# ---- Report ----
print("Number of nullcline crossings (steady states): %d" % len(steady_states))
for k, (xs, ys) in enumerate(steady_states, 1):
    fx = gX0 + gX1 / (1.0 + (ys / Yth) ** nY) - kX * xs
    fy = gY0 + gY1 / (1.0 + (xs / Xth) ** nX) - kY * ys
    print("Steady state %d: X = %.6f, Y = %.6f  (residual fX = %.3e, fY = %.3e)"
          % (k, xs, ys, fx, fy))

# ---- Check: exactly three crossings ----
three_crossings = (len(steady_states) == 3)
print("Exactly three crossings found:", three_crossings)
# Three crossings confirm bistability: because the two monotone (decreasing)
# nullclines meet an odd number of times, the outer two intersections are the
# two stable states and the middle one is the unstable saddle separating them.

# ---- Phase-plane plot ----
plt.figure(figsize=(8, 6))
plt.plot(X_of_Y, Y_sweep, 'b-', label='X-nullcline (dX/dt=0)')
plt.plot(X_sweep, Y_of_X, 'r-', label='Y-nullcline (dY/dt=0)')
for k, (xs, ys) in enumerate(steady_states, 1):
    plt.plot(xs, ys, 'ko', markersize=9)
    plt.annotate("SS%d (%.0f, %.0f)" % (k, xs, ys), (xs, ys),
                 textcoords="offset points", xytext=(8, 8))
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Toggle switch nullclines and steady states")
plt.legend()
plt.xlim(0, 600)
plt.ylim(0, 600)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.4.1_s5.png")
print("Figure saved.")
