import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# ----------------------------------------------------------------------
# Target: stationary distribution of the two-state bursting-gene promoter
# State 0 = inactive, State 1 = active.  We are given ONLY p_ss.
# ----------------------------------------------------------------------
p_ss = np.array([0.25, 0.75])          # target P(x) = (P(inactive), P(active))

def P(x):
    return p_ss[x]                      # target probability of state x

# ----------------------------------------------------------------------
# Metropolis-Hastings with a symmetric two-state proposal.
# The proposal always offers the OTHER state, so the proposal is symmetric
# and the acceptance ratio reduces to a = min(1, P(x')/P(x)).
# ----------------------------------------------------------------------
rng = np.random.default_rng(1)         # seed 1
n_steps = int(1e4)
x = 1                                  # start active (state 1)
chain = np.empty(n_steps + 1, dtype=int)
chain[0] = x

for t in range(n_steps):
    x_prop = 1 - x                     # symmetric proposal: always the other state
    a = min(1.0, P(x_prop) / P(x))     # acceptance probability
    if rng.random() < a:               # accept with probability a
        x = x_prop
    # else: reject and keep x
    chain[t + 1] = x                   # record the (possibly unchanged) state

# ----------------------------------------------------------------------
# Result 1: sampled state fractions (should approach p_ss).
# ----------------------------------------------------------------------
frac0 = np.mean(chain == 0)
frac1 = np.mean(chain == 1)
print(f"Sampled fraction inactive (state 0): {frac0:.4f}")
print(f"Sampled fraction active   (state 1): {frac1:.4f}")
print(f"Target stationary distribution p_ss: [{p_ss[0]:.4f}, {p_ss[1]:.4f}]")

# ----------------------------------------------------------------------
# Result 2: recover the sampler's transition matrix by counting the
# observed consecutive transitions in the realized chain.
# ----------------------------------------------------------------------
counts = np.zeros((2, 2))
for i in range(n_steps):
    counts[chain[i], chain[i + 1]] += 1
T_recovered = counts / counts.sum(axis=1, keepdims=True)   # row-normalize
print("Recovered MH transition matrix (rows = from-state):")
print(f"  T[0,0]={T_recovered[0,0]:.4f}  T[0,1]={T_recovered[0,1]:.4f}")
print(f"  T[1,0]={T_recovered[1,0]:.4f}  T[1,1]={T_recovered[1,1]:.4f}")

# Analytic MH transition matrix for comparison:
# from 0 propose 1: ratio P(1)/P(0)=3 -> accept 1
# from 1 propose 0: ratio P(0)/P(1)=1/3 -> accept 1/3
T_mh = np.array([[0.0, 1.0],
                 [1.0 / 3.0, 2.0 / 3.0]])
print("Analytic MH transition matrix:")
print(f"  T[0,0]={T_mh[0,0]:.4f}  T[0,1]={T_mh[0,1]:.4f}")
print(f"  T[1,0]={T_mh[1,0]:.4f}  T[1,1]={T_mh[1,1]:.4f}")

# ----------------------------------------------------------------------
# Separate check: a DIFFERENT chain (an example "original" 8B.1 promoter
# dynamics) that shares the same stationary distribution.
# Its off-diagonal rates satisfy the balance p_ss[0]*a = p_ss[1]*b.
# ----------------------------------------------------------------------
T_orig = np.array([[0.5, 0.5],
                   [1.0 / 6.0, 5.0 / 6.0]])   # 0.25*0.5 = 0.75*(1/6): same p_ss

def stationary(T):
    # left eigenvector of T with eigenvalue 1, normalized to sum 1
    vals, vecs = np.linalg.eig(T.T)
    v = np.real(vecs[:, np.argmin(np.abs(vals - 1.0))])
    return v / v.sum()

pi_mh = stationary(T_mh)
pi_orig = stationary(T_orig)
print("Stationary distribution of analytic MH chain:")
print(f"  [{pi_mh[0]:.4f}, {pi_mh[1]:.4f}]")
print("Stationary distribution of original 8B.1-style chain:")
print(f"  [{pi_orig[0]:.4f}, {pi_orig[1]:.4f}]")
print("Max |T_mh - T_orig| (transition matrices differ):")
print(f"  {np.max(np.abs(T_mh - T_orig)):.4f}")
print("Max |pi_mh - pi_orig| (same stationary distribution):")
print(f"  {np.max(np.abs(pi_mh - pi_orig)):.6f}")

# One-sentence explanation of why this check confirms the result:
print("Why the check confirms it: the MH chain and the original chain have "
      "different transition matrices yet identical stationary distributions, "
      "proving the sampler targets p_ss without needing the original dynamics, "
      "since many distinct chains can share one stationary distribution.")

# ----------------------------------------------------------------------
# Plot: sampled fractions vs target.
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4))
xpos = np.array([0, 1])
w = 0.35
ax.bar(xpos - w / 2, [frac0, frac1], width=w, label="MH sampled", color="steelblue")
ax.bar(xpos + w / 2, p_ss, width=w, label="target p_ss", color="darkorange")
ax.set_xticks(xpos)
ax.set_xticklabels(["inactive (0)", "active (1)"])
ax.set_ylabel("probability")
ax.set_title("Metropolis-Hastings reproduces the target stationary distribution")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.3.1_s1.png")
