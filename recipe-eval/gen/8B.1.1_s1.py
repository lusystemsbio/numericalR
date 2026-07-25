import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# Transition matrix: rows = current state (0=active, 1=inactive)
# T[0] = probabilities of next state given currently active
# T[1] = probabilities of next state given currently inactive
T = np.array([[0.7, 0.3],
              [0.1, 0.9]])

# --- Direct simulation of a two-state Markov chain ---
def simulate_promoter(T, start, n_steps, seed):
    rng = np.random.default_rng(seed)          # reproducible RNG
    states = np.empty(n_steps + 1, dtype=int)  # store trajectory (include start)
    states[0] = start                          # initial state
    for t in range(n_steps):
        u = rng.random()                       # draw a uniform in [0, 1)
        row = T[states[t]]                     # current row of T
        # pick next state: if u < P(stay in state 0-index 0) go to 0, else 1
        states[t + 1] = 0 if u < row[0] else 1
    return states

n_steps = int(1e4)
states = simulate_promoter(T, start=0, n_steps=n_steps, seed=1)  # start active (0)

# --- Running fractions of active / inactive over time ---
active_indicator = (states == 0).astype(float)   # 1 when active, 0 when inactive
cum_active = np.cumsum(active_indicator)          # cumulative active count
denom = np.arange(1, len(states) + 1)             # number of steps counted so far
running_active = cum_active / denom               # running active fraction
running_inactive = 1.0 - running_active           # running inactive fraction

# --- Numerical results ---
final_active = running_active[-1]
final_inactive = running_inactive[-1]
print(f"Final running active fraction:   {final_active:.4f}")
print(f"Final running inactive fraction: {final_inactive:.4f}")

# Analytic stationary distribution for reference: pi T = pi
# For this T, pi_active = 0.1/(0.3+0.1) = 0.25, pi_inactive = 0.75
pi_active = T[1, 0] / (T[0, 1] + T[1, 0])
pi_inactive = 1.0 - pi_active
print(f"Theoretical stationary active fraction:   {pi_active:.4f}")
print(f"Theoretical stationary inactive fraction: {pi_inactive:.4f}")

# --- Burst / spell check: count consecutive-run lengths in each state ---
run_lengths_active = []
run_lengths_inactive = []
cur = states[0]
count = 1
for s in states[1:]:
    if s == cur:
        count += 1
    else:
        (run_lengths_active if cur == 0 else run_lengths_inactive).append(count)
        cur = s
        count = 1
(run_lengths_active if cur == 0 else run_lengths_inactive).append(count)

mean_active_spell = np.mean(run_lengths_active)
mean_inactive_spell = np.mean(run_lengths_inactive)
print(f"Mean active spell length (steps):   {mean_active_spell:.4f}")
print(f"Mean inactive spell length (steps): {mean_inactive_spell:.4f}")
print(f"Theoretical mean active spell (1/0.3):   {1/0.3:.4f}")
print(f"Theoretical mean inactive spell (1/0.1): {1/0.1:.4f}")

# --- Plot ---
first = 200  # show promoter state over the first steps
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7))

ax1.step(np.arange(first + 1), states[:first + 1], where="post", color="tab:blue")
ax1.set_yticks([0, 1])
ax1.set_yticklabels(["active", "inactive"])
ax1.set_xlabel("step")
ax1.set_ylabel("promoter state")
ax1.set_title(f"Promoter state over first {first} steps (bursty: short active, long silent)")

ax2.plot(denom, running_active, label="running active fraction", color="tab:green")
ax2.plot(denom, running_inactive, label="running inactive fraction", color="tab:red")
ax2.axhline(0.25, ls="--", color="tab:green", alpha=0.6, label="0.25")
ax2.axhline(0.75, ls="--", color="tab:red", alpha=0.6, label="0.75")
ax2.set_xlabel("step")
ax2.set_ylabel("running fraction")
ax2.set_title("Running active/inactive fractions converging to 25% / 75%")
ax2.legend(loc="center right")

plt.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.1.1_s1.png")

# One-sentence explanation of why the check confirms the result:
print("Check explanation: The trajectory showing short active spells (~1/0.3 steps) "
      "and long inactive spells (~1/0.1 steps) with running fractions settling near "
      "0.25/0.75 matches the analytic stationary distribution pi=(0.25,0.75), "
      "confirming the simulation samples the correct long-run behavior.")
