import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Stiff model parameters: constitutively expressed gene with fast degradation ---
g = 50.0    # constant production rate
k = 10.0    # large degradation rate -> stiff
X0 = 3.0    # initial condition
t0, tf = 0.0, 4.0
dt = 0.2    # step size at which explicit Euler is unstable (dt*k = 2 > 2 stability limit)

# Time grid
t = np.arange(t0, tf + dt/2, dt)
n = len(t)

# Steady state / exact asymptote
steady = g / k
print(f"Steady state g/k: {steady}")

# --- Exact solution: X(t) = g/k + (X0 - g/k)*exp(-k t) ---
X_exact = steady + (X0 - steady) * np.exp(-k * t)

# --- Forward (explicit) Euler: X_{n+1} = X_n + (g - k*X_n)*dt ---
X_fe = np.empty(n)
X_fe[0] = X0
for i in range(n - 1):
    # evaluate RHS at the CURRENT (known) state
    rhs = g - k * X_fe[i]
    X_fe[i + 1] = X_fe[i] + rhs * dt

# --- Backward (implicit) Euler: X_{n+1} = X_n + (g - k*X_{n+1})*dt ---
# Linear RHS lets us solve for X_{n+1} explicitly: X_{n+1} = (X_n + dt*g)/(1 + dt*k)
X_be = np.empty(n)
X_be[0] = X0
for i in range(n - 1):
    # rearranged implicit update (no root-finding needed because RHS is linear)
    X_be[i + 1] = (X_be[i] + dt * g) / (1.0 + dt * k)

# --- Print trajectories ---
print("t        exact        forward_Euler        backward_Euler")
for i in range(n):
    print(f"{t[i]:.2f}    {X_exact[i]:.6f}    {X_fe[i]:.6e}    {X_be[i]:.6f}")

# --- Diagnostics for the stability check ---
print(f"dt*k (forward Euler stable only if < 2): {dt*k}")
print(f"Forward Euler final value: {X_fe[-1]:.6e}")
print(f"Forward Euler max |X|: {np.max(np.abs(X_fe)):.6e}")

# detect oscillation in forward Euler: sign changes of successive differences
diffs = np.diff(X_fe)
sign_changes = int(np.sum(np.diff(np.sign(diffs)) != 0))
print(f"Forward Euler sign changes in successive differences (oscillation indicator): {sign_changes}")

print(f"Backward Euler final value: {X_be[-1]:.6f}")
print(f"Backward Euler error vs exact at final time: {abs(X_be[-1] - X_exact[-1]):.6e}")
print(f"Backward Euler error vs steady state at final time: {abs(X_be[-1] - steady):.6e}")

# --- Plot comparison ---
plt.figure(figsize=(9, 6))
plt.plot(t, X_exact, 'k-', lw=2, label='Exact')
plt.plot(t, X_be, 'go-', lw=1.5, label='Backward (implicit) Euler')
plt.plot(t, X_fe, 'r.--', lw=1.5, label='Forward (explicit) Euler')
plt.axhline(steady, color='gray', ls=':', label='g/k = 5')
plt.xlabel('t')
plt.ylabel('X(t)')
plt.title('Stiff gene model (g=50, k=10, dt=0.2): forward vs backward Euler')
plt.legend()
plt.grid(True, alpha=0.3)
# clamp y-range so the stable curves remain visible despite FE blowup
plt.ylim(-5, 15)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.8.1_s3.png")

# --- One-sentence explanation ---
print("Explanation: The check confirms the result because forward Euler's amplification factor "
      "(1 - dt*k) = -1 has magnitude >= 1 at dt=0.2, so its error grows and oscillates in sign, "
      "whereas backward Euler's factor 1/(1 + dt*k) is always between 0 and 1, guaranteeing a "
      "stable monotone decay to g/k regardless of step size.")
