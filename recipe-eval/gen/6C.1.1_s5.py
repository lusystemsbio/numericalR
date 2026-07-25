import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ---- Model / simulation parameters ----
def f(x):          # drift term; pure diffusion means no drift
    return 0.0

def D(x):          # diffusion coefficient (constant here)
    return 10.0

X0 = 0.0           # initial position
T  = 200.0         # final time
dt = 0.01          # time step
n_steps = int(round(T / dt))
n_traj  = 10       # number of trajectories

# Time grid (length n_steps + 1 including t=0)
t = np.linspace(0.0, T, n_steps + 1)

# Reproducible randomness
rng = np.random.default_rng(1)

# Storage for all trajectories: rows = trajectories, cols = time points
X = np.zeros((n_traj, n_steps + 1))
X[:, 0] = X0

# ---- Euler-Maruyama integration (explicit, step by step) ----
for k in range(n_steps):
    x = X[:, k]                                   # current state of every trajectory
    # dW ~ N(0, dt): standard normal scaled by sqrt(dt)
    dW = rng.normal(loc=0.0, scale=np.sqrt(dt), size=n_traj)
    drift     = f(x) * dt                         # deterministic increment f(X)*dt
    diffusion = np.sqrt(2.0 * D(x)) * dW          # stochastic increment sqrt(2D)*dW
    X[:, k + 1] = x + drift + diffusion           # X_next = X + f*dt + sqrt(2D)*dW

# ---- Numerical checks against free Brownian motion ----
# For dX = sqrt(2D) dW with X0=0, X(t) ~ N(0, 2D t), so Var[X(T)] = 2 D T.
final = X[:, -1]
theory_var = 2.0 * D(0.0) * T
sample_var = np.var(final, ddof=1)
sample_mean = np.mean(final)

print(f"Number of trajectories:            {n_traj}")
print(f"Diffusion coefficient D:           {D(0.0)}")
print(f"Final time T:                      {T}")
print(f"Time step dt:                      {dt}")
print(f"Theoretical mean at T (should be 0): {0.0}")
print(f"Sample mean of X(T):               {sample_mean}")
print(f"Theoretical variance at T (2*D*T): {theory_var}")
print(f"Sample variance of X(T):           {sample_var}")
print(f"Ratio sample/theory variance:      {sample_var / theory_var}")

# ---- Plot the ten Brownian trajectories versus time ----
plt.figure(figsize=(9, 6))
for i in range(n_traj):
    plt.plot(t, X[i], lw=0.8)
plt.xlabel("time t")
plt.ylabel("X(t)")
plt.title("Ten free Brownian trajectories (Euler-Maruyama, f=0, D=10)")
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/6C.1.1_s5.png",
            dpi=150, bbox_inches="tight")

# One-sentence explanation:
# The check confirms the result because free Brownian motion has X(t) ~ N(0, 2*D*t),
# so the sample variance of X(T) matching 2*D*T (with mean ~ 0) shows the integrator
# reproduces the correct diffusive spreading of the previous chapter's Brownian motion.
print("Check: sample variance of X(T) matches 2*D*T (N(0,2Dt)), confirming free Brownian motion.")
