import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Toggle-switch parameters (genes X and Y mutually repress)
# ---------------------------------------------------------------
gX0, gX1, Yth, nY, kX = 5.0, 50.0, 100.0, 4.0, 0.1
gY0, gY1, Xth, nX, kY = 4.0, 40.0, 150.0, 4.0, 0.12

# ---------------------------------------------------------------
# Closed-form nullclines by separation of variables.
# fX(X,Y)=0  =>  gX0 + gX1/(1+(Y/Yth)^nY) - kX*X = 0
#            =>  X = (gX0 + gX1/(1+(Y/Yth)^nY)) / kX     (X explicit in Y)
# fY(X,Y)=0  =>  Y = (gY0 + gY1/(1+(X/Xth)^nX)) / kY     (Y explicit in X)
# ---------------------------------------------------------------
def X_nullcline(Y):          # X as a function of swept Y
    return (gX0 + gX1 / (1.0 + (Y / Yth) ** nY)) / kX

def Y_nullcline(X):          # Y as a function of swept X
    return (gY0 + gY1 / (1.0 + (X / Xth) ** nX)) / kY

# ---------------------------------------------------------------
# Sweep the free variable to trace each nullcline for plotting.
# ---------------------------------------------------------------
Y_sweep = np.linspace(0.0, 800.0, 2000)     # sweep Y -> get X on X-nullcline
X_on_Xnull = X_nullcline(Y_sweep)

X_sweep = np.linspace(0.0, 800.0, 2000)     # sweep X -> get Y on Y-nullcline
Y_on_Ynull = Y_nullcline(X_sweep)

# ---------------------------------------------------------------
# Find intersections (steady states) explicitly.
# At a crossing both hold, so substitute Y_nullcline(X) into the
# X-nullcline relation and look for roots of G(X) = X_nullcline(Y_nullcline(X)) - X.
# ---------------------------------------------------------------
def G(X):
    return X_nullcline(Y_nullcline(X)) - X

# Scan X, detect sign changes, refine each bracket by bisection.
Xscan = np.linspace(1e-6, 800.0, 4000)
Gvals = G(Xscan)
roots = []
for i in range(len(Xscan) - 1):
    if Gvals[i] == 0.0:
        roots.append(Xscan[i])
    elif Gvals[i] * Gvals[i + 1] < 0.0:          # sign change => a root inside
        a, b = Xscan[i], Xscan[i + 1]
        fa = Gvals[i]
        for _ in range(100):                     # bisection refinement
            m = 0.5 * (a + b)
            fm = G(m)
            if fa * fm <= 0.0:
                b = m
            else:
                a, fa = m, fm
        roots.append(0.5 * (a + b))

# Each root X* gives the matching Y* from the Y-nullcline.
steady_states = [(Xs, Y_nullcline(Xs)) for Xs in roots]

# ---------------------------------------------------------------
# Report numerical results.
# ---------------------------------------------------------------
print(f"Number of nullcline crossings (steady states): {len(steady_states)}")
for idx, (Xs, Ys) in enumerate(steady_states, 1):
    print(f"Steady state {idx}: X = {Xs:.6f}")
    print(f"Steady state {idx}: Y = {Ys:.6f}")
    # residuals confirm each point satisfies both dX/dt=0 and dY/dt=0
    resX = gX0 + gX1 / (1.0 + (Ys / Yth) ** nY) - kX * Xs
    resY = gY0 + gY1 / (1.0 + (Xs / Xth) ** nX) - kY * Ys
    print(f"Steady state {idx}: residual dX/dt = {resX:.3e}")
    print(f"Steady state {idx}: residual dY/dt = {resY:.3e}")

three_crossings = (len(steady_states) == 3)
print(f"Exactly three crossings found: {three_crossings}")

# ---------------------------------------------------------------
# Phase-plane plot of both nullclines and their crossings.
# ---------------------------------------------------------------
plt.figure(figsize=(7, 6))
plt.plot(X_on_Xnull, Y_sweep, 'b-', label='X-nullcline (dX/dt=0)')
plt.plot(X_sweep, Y_on_Ynull, 'r-', label='Y-nullcline (dY/dt=0)')
for idx, (Xs, Ys) in enumerate(steady_states, 1):
    plt.plot(Xs, Ys, 'ko', markersize=9)
    plt.annotate(f"SS{idx}", (Xs, Ys), textcoords="offset points", xytext=(8, 8))
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Toggle-switch nullclines and steady states')
plt.legend()
plt.xlim(0, 800)
plt.ylim(0, 800)
plt.grid(True, alpha=0.3)
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/3A.4.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Explanation: Each crossing is a point where both dX/dt=0 and dY/dt=0 "
      "simultaneously, i.e. a fixed point, so finding exactly three crossings "
      "confirms the toggle switch's bistable structure of two stable states "
      "flanking one unstable saddle seen in simulation.")
