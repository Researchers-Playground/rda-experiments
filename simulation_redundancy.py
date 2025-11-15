import numpy as np
import matplotlib.pyplot as plt
from math import exp, factorial

# Parameters
N = 5000
m = 100
lam = N / m

def poisson_cdf(k, lam):
    """Cumulative P(Y <= k) for Poisson(lam)."""
    return sum((lam**j / factorial(j)) for j in range(k+1)) * exp(-lam)

def prob_all_ge_X(X, lam, m):
    """Probability that all m columns have >= X nodes."""
    p_less = poisson_cdf(X-1, lam)
    return (1 - p_less)**m

# Compute probability for various X
X_values = np.arange(10, 41)  # adjust range for visibility
probs = [prob_all_ge_X(X, lam, m) for X in X_values]

# Plot
plt.figure(figsize=(7, 4))
plt.plot(X_values, probs, lw=2, marker='o', color='#1f77b4')
plt.axhline(0.99, color='red', linestyle='--', label='99% threshold')
plt.xlabel("Minimum honest nodes per column (X)")
plt.ylabel("Probability all 100 columns satisfy ≥ X")
plt.title("Redundancy Guarantee with N = 5000, m = 100 (λ = 50)")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig("redundancy_5000_nodes.png", dpi=300)
plt.show()