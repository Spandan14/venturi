import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
import sys
import os
import re
from matplotlib.widgets import Button

class PressureViewer:
    def __init__(self, directory):
        self.directory = directory
        self.current_idx = 0
        
        # Find all pressure matrix files
        pressure_files = []
        for f in os.listdir(directory):
            if f.startswith('pressure_matrix_') and f.endswith('.mtx'):
                match = re.search(r'pressure_matrix_([0-9.]+)\.mtx', f)
                if match:
                    timestamp = float(match.group(1))
                    pressure_files.append((timestamp, f))
        
        # Sort by timestamp
        pressure_files.sort(key=lambda x: x[0])
        self.timestamps = [t for t, _ in pressure_files]
        self.pressure_files = [f for _, f in pressure_files]
        
        # Find corresponding divergence files
        self.div_files = []
        for timestamp in self.timestamps:
            div_file = f'divergence_{timestamp}.vec'
            if os.path.exists(os.path.join(directory, div_file)):
                self.div_files.append(div_file)
            else:
                print(f"Warning: No divergence file found for timestamp {timestamp}")
                self.div_files.append(None)
        
        if len(self.pressure_files) == 0:
            print(f"No pressure matrix files found in {directory}")
            sys.exit(1)
        
        print(f"Found {len(self.pressure_files)} timesteps")
        print(f"Time range: {self.timestamps[0]} to {self.timestamps[-1]}")
        
        # Create figure
        self.fig, self.axes = plt.subplots(1, 3, figsize=(16, 5))
        self.fig.subplots_adjust(bottom=0.15)
        self.colorbar = None
        
        # Add navigation buttons
        ax_prev = plt.axes([0.3, 0.02, 0.1, 0.05])
        ax_next = plt.axes([0.6, 0.02, 0.1, 0.05])
        self.btn_prev = Button(ax_prev, 'Previous')
        self.btn_next = Button(ax_next, 'Next')
        self.btn_prev.on_clicked(self.prev_frame)
        self.btn_next.on_clicked(self.next_frame)
        
        # Display first frame
        self.update_plot()
        
        # Keyboard navigation (after initial plot)
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        plt.show()
    
    def load_matrix(self, filename):
        filepath = os.path.join(self.directory, filename)
        data = []
        row_ind = []
        col_ind = []
        max_row = 0
        max_col = 0
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        row = int(parts[0])
                        col = int(parts[1])
                        val = float(parts[2])
                        if row >= 0 and col >= 0:
                            row_ind.append(row)
                            col_ind.append(col)
                            data.append(val)
                            max_row = max(max_row, row)
                            max_col = max(max_col, col)
                    except (ValueError, IndexError):
                        continue
        
        return csr_matrix((data, (row_ind, col_ind)), shape=(max_row+1, max_col+1))
    
    def load_divergence(self, filename):
        if filename is None:
            return None
        filepath = os.path.join(self.directory, filename)
        if os.path.exists(filepath):
            return np.loadtxt(filepath)
        return None
    
    def update_plot(self):
        # Remove old colorbar if it exists
        if hasattr(self, 'colorbar') and self.colorbar is not None:
            self.colorbar.remove()
            self.colorbar = None
        
        # Clear previous plots
        for ax in self.axes:
            ax.clear()
        
        timestamp = self.timestamps[self.current_idx]
        A = self.load_matrix(self.pressure_files[self.current_idx])
        divergence = self.load_divergence(self.div_files[self.current_idx])
        
        # 1. Spy plot (structure)
        self.axes[0].spy(A, markersize=1)
        self.axes[0].set_title(f'Matrix Structure\nt={timestamp:.6f}')
        self.axes[0].set_xlabel('Column')
        self.axes[0].set_ylabel('Row')
        
        # 2. Full matrix heatmap
        im = self.axes[1].imshow(A.toarray(), cmap='RdBu_r', aspect='auto')
        self.axes[1].set_title('Matrix Values')
        self.axes[1].set_xlabel('Column')
        self.axes[1].set_ylabel('Row')
        self.colorbar = plt.colorbar(im, ax=self.axes[1])
        
        # 3. Divergence vector
        if divergence is not None:
            self.axes[2].plot(divergence, 'o-', markersize=2)
            self.axes[2].set_title('Divergence Vector')
            self.axes[2].set_xlabel('Index')
            self.axes[2].set_ylabel('Value')
            self.axes[2].grid(True, alpha=0.3)
            
            # Print stats
            print(f"\nFrame {self.current_idx+1}/{len(self.timestamps)} | t={timestamp:.6f}")
            print(f"  Matrix: {A.shape[0]}x{A.shape[1]}, {A.nnz} non-zeros, "
                  f"{100 * (1 - A.nnz / (A.shape[0] * A.shape[1])):.2f}% sparse")
            print(f"  Divergence range: [{divergence.min():.6e}, {divergence.max():.6e}]")
        else:
            self.axes[2].text(0.5, 0.5, 'No divergence data', 
                            ha='center', va='center', transform=self.axes[2].transAxes)
        
        self.fig.suptitle(f'Timestep {self.current_idx+1}/{len(self.timestamps)} | '
                         f'Use arrow keys or buttons to navigate', fontsize=12)
        self.fig.canvas.draw()
    
    def next_frame(self, event=None):
        if self.current_idx < len(self.timestamps) - 1:
            self.current_idx += 1
            self.update_plot()
    
    def prev_frame(self, event=None):
        if self.current_idx > 0:
            self.current_idx -= 1
            self.update_plot()
    
    def on_key(self, event):
        print(f"Key pressed: {event.key}")  # Debug output
        if event.key == 'right' or event.key == 'n':
            self.next_frame()
        elif event.key == 'left' or event.key == 'p':
            self.prev_frame()
        elif event.key == 'q':
            plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python visualize_pressure.py <directory>")
        print("Example: python visualize_pressure.py ./output")
        print("\nExpected file format:")
        print("  pressure_matrix_<timestamp>.mtx")
        print("  divergence_<timestamp>.txt")
        sys.exit(1)
    
    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory")
        sys.exit(1)
    
    viewer = PressureViewer(directory)
