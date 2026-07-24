import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- Model parameters (toggle switch) ----------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1   # dX/dt = gX0 + gX1/(1+(Y/Yth)^nY) - kX*X
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12  # dY/dt = gY0 + gY1/(1+(X/Xth)^nX) - kY*Y

# ---------------- Closed-form nullclines by separation of variables ----------------
# X-nullcline: set fX(X,Y)=0 and solve for X (X appears only linearly via -kX*X):
#   0 = gX0 + gX1/(1+(Y/Yth)^nY) - kX*X  ->  X = [gX0 + gX1/(1+(Y/Yth)^nY)] / kX
def X_of_Y(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# Y-nullcline: set fY(X,Y)=0 and solve for Y (Y appears only linearly via -kY*Y):
#   0 = gY0 + gY1/(1+(X/Xth)^nX) - kY*Y  ->  Y = [gY0 + gY1/(1+(X/Xth)^nX)] / kY
def Y_of_X(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# Sweep the "other" variable to trace each nullcline as a curve in the (X,Y) plane
Y_sweep = np.linspace(0, 600, 2000)   # for X-nullcline we sweep Y
X_on_Xnull = X_of_Y(Y_sweep)          # matching X values

X_sweep = np.linspace(0, 600, 2000)   # for Y-nullcline we sweep X
Y_on_Ynull = Y_of_X(X_sweep)          # matching Y values

# ---------------- Find the intersections (steady states) explicitly ----------------
# A steady state satisfies BOTH nullclines. Compose them: given X, the Y-nullcline
# fixes Y = Y_of_X(X); feeding that into the X-nullcline should return the same X.
# So define g(X) = X_of_Y(Y_of_X(X)) - X and hunt for its sign changes (roots).
def g(X):
    return X_of_Y(Y_of_X(X)) - X

X_scan = np.linspace(1e-6, 600, 6000)
gvals = g(X_scan)

roots = []
for i in range(len(X_scan) - 1):
    if gvals[i] == 0.0:
        roots.append(X_scan[i])
    elif gvals[i] * gvals[i + 1] < 0.0:
        # bisection to refine the bracketed root
        a, b = X_scan[i], X_scan[i + 1]
        for _ in range(80):
            m = 0.5 * (a + b)
            if g(a) * g(m) <= 0.0:
                b = m
            else:
                a = m
        roots.append(0.5 * (a + b))

# Deduplicate near-identical roots
roots_clean = []
for r in roots:
    if not any(abs(r - s) < 1e-3 for s in roots_clean):
        roots_clean.append(r)

steady_states = [(r, Y_of_X(r)) for r in roots_clean]

# ---------------- Report numerical results ----------------
print(f"Number of nullcline crossings (steady states) found: {len(steady_states)}")
for k, (xs, ys) in enumerate(steady_states, 1):
    fX = gX0 + gX1 / (1.0 + (ys / Yth) ** nY) - kX * xs
    fY = gY0 + gY1 / (1.0 + (xs / Xth) ** nX) - kY * ys
    print(f"Steady state {k}: X = {xs:.6f}, Y = {ys:.6f}  (residual fX = {fX:.3e}, fY = {fY:.3e})")

print(f"Three-crossing check passed: {len(steady_states) == 3}")

# ---------------- Phase-plane plot ----------------
plt.figure(figsize=(7, 6))
plt.plot(X_on_Xnull, Y_sweep, 'b-', label='X-nullcline (dX/dt=0)')
plt.plot(X_sweep, Y_on_Ynull, 'r-', label='Y-nullcline (dY/dt=0)')
for xs, ys in steady_states:
    plt.plot(xs, ys, 'ko', markersize=9)
    plt.annotate(f"({xs:.0f},{ys:.0f})", (xs, ys),
                 textcoords="offset points", xytext=(8, 8))
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Toggle switch nullclines and steady states")
plt.legend()
plt.xlim(0, 600)
plt.ylim(0, 600)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.4.1_s3.png")

# One-sentence explanation of why the check confirms the result:
# The nullclines are the loci where each variable stops changing, so every point
# where they cross is a simultaneous root of both ODEs (a steady state); finding
# exactly three crossings confirms the bistable toggle-switch structure of two
# stable states flanking one unstable state, matching what simulation shows.
print("Explanation: crossings of the two nullclines are exactly the points where "
      "both dX/dt=0 and dY/dt=0 simultaneously, so three crossings confirm three "
      "steady states (two stable, one unstable) as seen by simulation.")
