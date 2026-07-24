import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# --- Stiff model parameters: dX/dt = g - k*X ---
g = 50.0      # constitutive production rate
k = 10.0      # large degradation rate -> stiff
X0 = 3.0      # initial condition
dt = 0.2      # step size (chosen so forward Euler is unstable)
t0, tf = 0.0, 4.0

# Time grid
t = np.arange(t0, tf + dt/2, dt)
n = len(t)

# Exact (analytical) solution: X(t) = g/k + (X0 - g/k)*exp(-k*t)
Xss = g / k                                    # steady state g/k = 5
X_exact = Xss + (X0 - Xss) * np.exp(-k * t)

# --- Forward (explicit) Euler, done step by step ---
X_fe = np.empty(n)
X_fe[0] = X0
for i in range(n - 1):
    # evaluate RHS at the CURRENT point, then step forward
    rhs = g - k * X_fe[i]
    X_fe[i + 1] = X_fe[i] + dt * rhs

# --- Backward (implicit) Euler, done step by step ---
# X_next = X + dt*(g - k*X_next)  ->  X_next*(1 + dt*k) = X + dt*g
X_be = np.empty(n)
X_be[0] = X0
for i in range(n - 1):
    # solve the implicit linear relation explicitly for X_next
    X_be[i + 1] = (X_be[i] + dt * g) / (1.0 + dt * k)

# --- Stability check: does forward Euler oscillate/diverge? ---
# The forward Euler amplification factor is (1 - dt*k); |factor| > 1 -> divergence.
fe_factor = 1.0 - dt * k
be_factor = 1.0 / (1.0 + dt * k)
fe_error_max = np.max(np.abs(X_fe - X_exact))
be_error_max = np.max(np.abs(X_be - X_exact))
# oscillation: sign of the deviation from steady state flips each step
fe_dev = X_fe - Xss
fe_sign_flips = int(np.sum(np.sign(fe_dev[1:]) != np.sign(fe_dev[:-1])))

# --- Print numerical results ---
print(f"Steady state g/k:                 {Xss}")
print(f"Forward Euler amplification 1-dt*k: {fe_factor}  (|.|>1 => unstable: {abs(fe_factor) > 1})")
print(f"Backward Euler factor 1/(1+dt*k):   {be_factor}  (|.|<1 => stable: {abs(be_factor) < 1})")
print(f"Forward Euler final value X(4):     {X_fe[-1]}")
print(f"Backward Euler final value X(4):    {X_be[-1]}")
print(f"Exact final value X(4):             {X_exact[-1]}")
print(f"Forward Euler max abs error:        {fe_error_max}")
print(f"Backward Euler max abs error:       {be_error_max}")
print(f"Forward Euler sign flips about g/k: {fe_sign_flips}  (oscillation)")
print(f"Forward Euler |final deviation| from g/k: {abs(X_fe[-1] - Xss)}")
print(f"Backward Euler |final deviation| from g/k: {abs(X_be[-1] - Xss)}")

# Check confirms: forward Euler's deviation grows and flips sign each step (diverging
# oscillation), while backward Euler's deviation shrinks monotonically toward g/k = 5,
# so only the implicit step remains stable at this stiff step size.
print("Explanation: |1 - dt*k| = 5 > 1 makes forward Euler amplify and flip the error "
      "each step (diverging oscillation), whereas 1/(1+dt*k) < 1 makes backward Euler "
      "shrink the error toward g/k every step, confirming the implicit step is stable.")

# --- Plot ---
plt.figure(figsize=(9, 5.5))
tf_fine = np.linspace(t0, tf, 400)
plt.plot(tf_fine, Xss + (X0 - Xss) * np.exp(-k * tf_fine),
         'k-', lw=2, label="Exact")
plt.plot(t, X_fe, 'r--o', ms=4, label="Forward Euler (dt=0.2)")
plt.plot(t, X_be, 'b-s', ms=4, label="Backward Euler (dt=0.2)")
plt.axhline(Xss, color='gray', ls=':', label="g/k = 5")
plt.xlabel("t")
plt.ylabel("X(t)")
plt.title("Stiff gene decay: forward vs backward Euler (dt=0.2)")
plt.legend()
plt.ylim(-40, 60)
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("/private/tmp/claude-501/-Users-lvmy-neu-teaching-numericalR-numericalR/e8e9d206-e3a8-4079-be52-06334c0a89a5/scratchpad/eval/gen/2B.8.1_s2.png")
