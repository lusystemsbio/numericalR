import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------
# Target stationary distribution of the two-state bursting-gene
# promoter (state 0 = active, state 1 = inactive).
# ---------------------------------------------------------------
P = np.array([0.25, 0.75])   # target P(x): p_ss = (0.25, 0.75)

# ---------------------------------------------------------------
# Metropolis-Hastings, implemented explicitly.
# Symmetric proposal: from any state always propose the OTHER state.
# Acceptance prob: a = min(1, P(x')/P(x)).
# ---------------------------------------------------------------
rng = np.random.default_rng(1)   # seed 1
n_steps = int(1e4)

x = 0                            # start active (state 0)
samples = np.empty(n_steps + 1, dtype=int)
samples[0] = x
for t in range(n_steps):
    x_prop = 1 - x                       # propose the other state (symmetric)
    a = min(1.0, P[x_prop] / P[x])       # acceptance ratio
    if rng.random() < a:                 # accept with prob a
        x = x_prop
    # else: stay in current state (implicit)
    samples[t + 1] = x

# ---------------------------------------------------------------
# Sampled state fractions
# ---------------------------------------------------------------
frac0 = np.mean(samples == 0)
frac1 = np.mean(samples == 1)
print(f"Target stationary distribution p_ss:      [{P[0]:.4f}, {P[1]:.4f}]")
print(f"Sampled fraction in state 0 (active):     {frac0:.4f}")
print(f"Sampled fraction in state 1 (inactive):   {frac1:.4f}")

# ---------------------------------------------------------------
# Recovered transition matrix from the empirical trajectory
# T_recovered[i, j] = fraction of transitions from i that go to j
# ---------------------------------------------------------------
counts = np.zeros((2, 2))
for i, j in zip(samples[:-1], samples[1:]):
    counts[i, j] += 1
T_recovered = counts / counts.sum(axis=1, keepdims=True)
print("Recovered transition matrix (from samples):")
print(f"  T[0->0]={T_recovered[0,0]:.4f}  T[0->1]={T_recovered[0,1]:.4f}")
print(f"  T[1->0]={T_recovered[1,0]:.4f}  T[1->1]={T_recovered[1,1]:.4f}")

# ---------------------------------------------------------------
# Analytic MH transition matrix (what the sampler is drawing from):
#   from 0: propose 1, a=min(1,0.75/0.25)=1        -> 0->1 = 1
#   from 1: propose 0, a=min(1,0.25/0.75)=1/3      -> 1->0 = 1/3, 1->1 = 2/3
# ---------------------------------------------------------------
T_MH = np.array([[0.0, 1.0],
                 [1.0/3.0, 2.0/3.0]])
print("Analytic MH transition matrix:")
print(f"  T[0->0]={T_MH[0,0]:.4f}  T[0->1]={T_MH[0,1]:.4f}")
print(f"  T[1->0]={T_MH[1,0]:.4f}  T[1->1]={T_MH[1,1]:.4f}")

# ---------------------------------------------------------------
# Separate check: an ORIGINAL bursting-gene chain (8B.1) with a
# DIFFERENT transition matrix but the SAME stationary distribution.
# Any chain with P(0->1)/P(1->0) = p_ss[1]/p_ss[0]*... i.e. ratio 3.
# ---------------------------------------------------------------
T_orig = np.array([[0.50, 0.50],
                   [1.0/6.0, 5.0/6.0]])   # different chain, still bursting-gene

def stationary(T):
    # left eigenvector with eigenvalue 1, normalized
    vals, vecs = np.linalg.eig(T.T)
    v = np.real(vecs[:, np.argmin(np.abs(vals - 1.0))])
    return v / v.sum()

pi_MH = stationary(T_MH)
pi_orig = stationary(T_orig)
print("Original chain (8B.1) transition matrix:")
print(f"  T[0->0]={T_orig[0,0]:.4f}  T[0->1]={T_orig[0,1]:.4f}")
print(f"  T[1->0]={T_orig[1,0]:.4f}  T[1->1]={T_orig[1,1]:.4f}")
print(f"Stationary dist of MH chain:        [{pi_MH[0]:.4f}, {pi_MH[1]:.4f}]")
print(f"Stationary dist of original chain:  [{pi_orig[0]:.4f}, {pi_orig[1]:.4f}]")
matrices_differ = not np.allclose(T_MH, T_orig)
stationaries_match = np.allclose(pi_MH, pi_orig)
print(f"Transition matrices differ:          {matrices_differ}")
print(f"Stationary distributions match:      {stationaries_match}")

# The two chains have different transition matrices yet the same stationary
# distribution, confirming that many distinct chains share one stationary
# distribution -- so MH's own chain need not equal the original to reproduce p_ss.

# ---------------------------------------------------------------
# Plot: sampled vs target
# ---------------------------------------------------------------
fig, ax = plt.subplots(figsize=(6, 4))
xpos = np.array([0, 1])
w = 0.35
ax.bar(xpos - w/2, [frac0, frac1], width=w, label="MH samples", color="steelblue")
ax.bar(xpos + w/2, P, width=w, label="target p_ss", color="orange")
ax.set_xticks(xpos)
ax.set_xticklabels(["active (0)", "inactive (1)"])
ax.set_ylabel("probability")
ax.set_title("Metropolis-Hastings reproduces the target stationary distribution")
ax.legend()
fig.tight_layout()
plt.savefig("/Users/lvmy/neu/teaching/numericalR/numericalR/recipe-eval/gen/8B.3.1_s3.png")
