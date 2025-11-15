import matplotlib.pyplot as plt
import numpy as np

# Data
schemes = ['RDA', 'CRDA (k=8)', 'CRDA (k=16)']
rep_factors = [50, 6.25, 3.125]
colors = ['#4472C4', '#ED7D31', '#A5A5A5']

# Plot
plt.figure(figsize=(6,4))
bars = plt.bar(schemes, rep_factors, color=colors, width=0.55)

# Add text labels
for bar in bars:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 1, f"{yval:.1f}×", 
             ha='center', va='bottom', fontsize=11, fontweight='medium')

# Labels & style
plt.ylabel("Replication Factor (× Block Size)")
plt.title("Reduction of Replication Factor via RLNC")
plt.ylim(0, 60)
plt.grid(axis='y', linestyle='--', alpha=0.4)
plt.tight_layout()

# Save and show
plt.show()