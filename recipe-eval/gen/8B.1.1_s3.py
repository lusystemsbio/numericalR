import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Two-state promoter Markov chain: state 0 = active, state 1 = inactive
# Rows of T give the transition probabilities FROM the current state.
T = np.array([[0.7, 0.3],   # from active:   stay active 0.7, leave to inactive 0.3
              [0.1, 0.9]])   # from inactive: to active 0.1, stay inactive 0.9

n_steps = int(1e4)
rng = np.random.default_rng(1)   # seed 1 for reproducibility

# --- Direct simulation: draw a uniform, pick next state from the current row ---
states = np.empty(n_steps, dtype=int)
states[0] = 0                     # start in the active state
for k in range(1, n_steps):
    row = T[states[k - 1]]        # transition probabilities from the current state
    u = rng.random()             # uniform draw in [0, 1)
    # cumulative comparison: next state is the first bin whose cumsum exceeds u
    states[k] = 0 if u < row[0] else 1

# --- Running fractions of time spent active / inactive ---
step_index = np.arange(1, n_steps + 1)
running_active = np.cumsum(states == 0) / step_index
running_inactive = np.cumsum(states == 1) / step_index

final_active = running_active[-1]
final_inactive = running_inactive[-1]

# --- Analytic stationary distribution for comparison ---
# Solve pi = pi @ T (left eigenvector for eigenvalue 1), normalized to sum 1.
vals, vecs = np.linalg.eig(T.T)
pi = np.real(vecs[:, np.argmin(np.abs(vals - 1))])
pi = pi / pi.sum()

print(f"Final running fraction ACTIVE:   {final_active:.4f}")
print(f"Final running fraction INACTIVE: {final_inactive:.4f}")
print(f"Analytic stationary fraction ACTIVE:   {pi[0]:.4f}")
print(f"Analytic stationary fraction INACTIVE: {pi[1]:.4f}")

# --- Burst check: measure lengths of consecutive active and inactive spells ---
change_points = np.where(np.diff(states) != 0)[0] + 1
segment_bounds = np.concatenate(([0], change_points, [n_steps]))
seg_lengths = np.diff(segment_bounds)
seg_states = states[segment_bounds[:-1]]
mean_active_spell = seg_lengths[seg_states == 0].mean()
mean_inactive_spell = seg_lengths[seg_states == 1].mean()

print(f"Mean ACTIVE spell length (steps):   {mean_active_spell:.4f}")
print(f"Mean INACTIVE spell length (steps): {mean_inactive_spell:.4f}")
# Explanation: because the active state leaves with prob 0.3 (short spells) and the
# inactive with prob 0.1 (long spells), the running fraction settling at ~25% active /
# ~75% inactive and matching the analytic stationary distribution confirms the chain
# reproduces the expected bursty steady state.

# --- Plots ---
window = 300   # first steps, so bursts are visible
fig, axes = plt.subplots(2, 1, figsize=(10, 7))

axes[0].step(np.arange(window), states[:window], where="post", color="tab:blue")
axes[0].set_yticks([0, 1])
axes[0].set_yticklabels(["active", "inactive"])
axes[0].set_xlabel("step")
axes[0].set_ylabel("promoter state")
axes[0].set_title(f"Promoter state over first {window} steps (bursty)")

axes[1].plot(step_index, running_active, label="running active fraction", color="tab:green")
axes[1].plot(step_index, running_inactive, label="running inactive fraction", color="tab:red")
axes[1].axhline(0.25, ls="--", color="gray", lw=1)
axes[1].axhline(0.75, ls="--", color="gray", lw=1)
axes[1].set_xlabel("step")
axes[1].set_ylabel("fraction")
axes[1].set_title("Running active / inactive fractions converge to 25% / 75%")
axes[1].legend()

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.1.1_s3.png")
