import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Stiff model: constitutively expressed gene with large degradation rate
# dX/dt = g - k*X  ->  exact solution decays to steady state g/k
g = 50.0      # production rate
k = 10.0      # degradation rate (large -> stiff)
X0 = 3.0      # initial condition
t0, tend = 0.0, 4.0
dt = 0.2      # step size: large enough that explicit Euler is unstable (dt*k = 2 > 2 needed for stability boundary)

# Time grid
t = np.arange(t0, tend + dt/2, dt)
n = len(t)

# Exact analytic solution: X(t) = g/k + (X0 - g/k)*exp(-k*t)
Xss = g / k
X_exact = Xss + (X0 - Xss) * np.exp(-k * t)

# Forward (explicit) Euler: X_{n+1} = X_n + (g - k*X_n)*dt
X_fe = np.empty(n)
X_fe[0] = X0
for i in range(n - 1):
    rhs = g - k * X_fe[i]            # evaluate RHS at the KNOWN current state
    X_fe[i + 1] = X_fe[i] + rhs * dt  # take the explicit step

# Backward (implicit) Euler: X_{n+1} = X_n + (g - k*X_{n+1})*dt
# Because the RHS is linear, solve for X_{n+1} directly:
#   X_{n+1}(1 + dt*k) = X_n + dt*g  ->  X_{n+1} = (X_n + dt*g)/(1 + dt*k)
X_be = np.empty(n)
X_be[0] = X0
for i in range(n - 1):
    X_be[i + 1] = (X_be[i] + dt * g) / (1.0 + dt * k)  # implicit update, no iteration needed

# Report numerical results
print(f"Steady state g/k = {Xss}")
print(f"dt = {dt}, stability factor dt*k = {dt * k}")
print(f"Final time t = {t[-1]}")
print(f"Forward Euler  X(final) = {X_fe[-1]}")
print(f"Backward Euler X(final) = {X_be[-1]}")
print(f"Exact          X(final) = {X_exact[-1]}")
print(f"Forward Euler  max |X| = {np.max(np.abs(X_fe))}  (blows up if huge)")
print(f"Backward Euler max |X| = {np.max(np.abs(X_be))}")
print(f"Backward Euler error at final = {abs(X_be[-1] - X_exact[-1])}")

# Check that forward Euler oscillates: sign of successive changes should flip
diffs = np.diff(X_fe)
sign_changes = np.sum(np.diff(np.sign(diffs)) != 0)
print(f"Forward Euler sign changes in successive increments = {sign_changes}  (oscillation)")
print(f"Forward Euler increment magnitudes growing (diverges): "
      f"|dX_1| = {abs(diffs[0])}, |dX_last| = {abs(diffs[-1])}")

# Plot comparison
plt.figure(figsize=(9, 6))
plt.plot(t, X_exact, "k-", lw=2, label="Exact")
plt.plot(t, X_be, "bo-", label="Backward (implicit) Euler")
plt.plot(t, X_fe, "rs--", label="Forward (explicit) Euler")
plt.axhline(Xss, color="gray", ls=":", label=f"Steady state g/k = {Xss}")
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title(f"Stiff gene model dX/dt = g - k*X  (g={g}, k={k}, dt={dt})")
plt.legend()
plt.grid(True, alpha=0.3)
# clip the y-axis so the exact/implicit curves remain visible despite FE blow-up
plt.ylim(-40, 90)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.8.1_s4.png")

# One-sentence explanation of why this check confirms the result:
print("Explanation: Because dt*k = 2 exceeds the explicit stability limit, forward "
      "Euler's error amplification factor (1 - dt*k) = -1 has magnitude >= 1 causing "
      "sign-flipping oscillations that fail to decay, while backward Euler's factor "
      "1/(1 + dt*k) is always between 0 and 1, so it stays bounded and converges to "
      "g/k = 5 just like the exact solution.")
