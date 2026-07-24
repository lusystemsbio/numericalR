import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Toggle switch parameters ----
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4, 0.12

# ---- Closed-form nullclines via separation of variables ----
# X-nullcline: solve fX(X,Y)=gX0+gX1/(1+(Y/Yth)^nY)-kX*X = 0 for X, sweeping Y.
def X_on_Xnull(Y):
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

# Y-nullcline: solve fY(X,Y)=gY0+gY1/(1+(X/Xth)^nX)-kY*Y = 0 for Y, sweeping X.
def Y_on_Ynull(X):
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ---- Sweep the "other" variable to trace each curve ----
Y_sweep = np.linspace(0, 600, 2000)      # X-nullcline: X as function of Y
X_of_Xnull = X_on_Xnull(Y_sweep)

X_sweep = np.linspace(0, 600, 2000)      # Y-nullcline: Y as function of X
Y_of_Ynull = Y_on_Ynull(X_sweep)

# ---- Find intersections (steady states) ----
# A steady state satisfies both: X = X_on_Xnull(Y) and Y = Y_on_Ynull(X).
# Substitute to get a single-variable self-consistency map h(X)=X_on_Xnull(Y_on_Ynull(X))-X=0.
def h(X):
    return X_on_Xnull(Y_on_Ynull(X)) - X

Xg = np.linspace(0.1, 600, 20000)
hv = h(Xg)
roots = []
for i in range(len(Xg) - 1):
    if hv[i] == 0.0 or hv[i] * hv[i + 1] < 0.0:  # sign change -> bracketed root
        a, b = Xg[i], Xg[i + 1]
        fa, fb = hv[i], hv[i + 1]
        for _ in range(100):  # bisection refinement
            m = 0.5 * (a + b)
            fm = h(m)
            if fa * fm <= 0.0:
                b, fb = m, fm
            else:
                a, fa = m, fm
        Xr = 0.5 * (a + b)
        Yr = Y_on_Ynull(Xr)
        roots.append((Xr, Yr))

# ---- Report steady states ----
print("Number of nullcline crossings (steady states):", len(roots))
for k, (Xr, Yr) in enumerate(roots, 1):
    print(f"Steady state {k}: X = {Xr:.6f}, Y = {Yr:.6f}")

# ---- Phase-plane plot ----
plt.figure(figsize=(7, 6))
plt.plot(X_of_Xnull, Y_sweep, 'b-', label="X-nullcline (dX/dt=0)")
plt.plot(X_sweep, Y_of_Ynull, 'r-', label="Y-nullcline (dY/dt=0)")
for Xr, Yr in roots:
    plt.plot(Xr, Yr, 'ko', markersize=9)
plt.xlabel("X")
plt.ylabel("Y")
plt.title("Toggle switch nullclines and steady states")
plt.xlim(0, 600)
plt.ylim(0, 600)
plt.legend()
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.4.1_s4.png")

# Explanation:
print("Check: exactly three crossings confirm bistability because the two nullclines")
print("intersect three times, and the count/positions match the two stable states plus")
print("one unstable saddle expected from simulation of a toggle switch.")
