import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import sys

# Get filepaths from command line arguments
if len(sys.argv) < 3:
    print("Usage: python visualize_pressure.py <divergence_file> <matrix_file>")
    print("Example: python visualize_pressure.py divergence.txt pressure_matrix.txt")
    sys.exit(1)

divergence_file = sys.argv[1]
matrix_file = sys.argv[2]

# Read the divergence vector
divergence = np.loadtxt(divergence_file)

# Read the sparse matrix
# Eigen outputs sparse matrices in a readable format
with open(matrix_file, 'r') as f:
    lines = f.readlines()

# Parse the matrix - Eigen format shows non-zero entries
data = []
row_ind = []
col_ind = []
max_row = 0
max_col = 0

for line in lines:
    parts = line.strip().split()
    if len(parts) >= 3:
        try:
            # Try to parse as triplet (row, col, value)
            row = int(parts[0])
            col = int(parts[1])
            val = float(parts[2])
            row_ind.append(row)
            col_ind.append(col)
            data.append(val)
            max_row = max(max_row, row)
            max_col = max(max_col, col)
        except ValueError:
            continue

# Create sparse matrix
A = csr_matrix((data, (row_ind, col_ind)), shape=(max_row+1, max_col+1))

# Create visualization
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# 1. Spy plot (structure)
axes[0].spy(A, markersize=1)
axes[0].set_title('Matrix Structure (Non-zero Pattern)')
axes[0].set_xlabel('Column')
axes[0].set_ylabel('Row')

# 2. Full matrix heatmap (dense version)
axes[1].imshow(A.toarray(), cmap='RdBu_r', aspect='auto')
axes[1].set_title('Matrix Values (Heatmap)')
axes[1].set_xlabel('Column')
axes[1].set_ylabel('Row')
plt.colorbar(axes[1].images[0], ax=axes[1])

# 3. Divergence vector
axes[2].plot(divergence, 'o-', markersize=2)
axes[2].set_title('Divergence Vector')
axes[2].set_xlabel('Index')
axes[2].set_ylabel('Value')
axes[2].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('pressure_system_visualization.png', dpi=150)
plt.show()

# Print matrix statistics
print(f"Matrix size: {A.shape[0]} x {A.shape[1]}")
print(f"Non-zeros: {A.nnz}")
print(f"Sparsity: {100 * (1 - A.nnz / (A.shape[0] * A.shape[1])):.2f}%")
print(f"Divergence range: [{divergence.min():.6f}, {divergence.max():.6f}]")
