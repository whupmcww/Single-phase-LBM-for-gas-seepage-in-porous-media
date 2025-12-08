import numpy as np
import matplotlib.pyplot as plt

size=1000
width=10
spacing=30

def generate_straight_channels(size, width, spacing):
    """Generate vertical straight channels and ensure fully contained within the range"""
    grid = np.ones((size, size))
    
    # Calculate the number of flow channels that can be accommodated
    channel_count = (size - width) // spacing
    start_pos = (size - (channel_count * spacing + width)) // 2  # Start position
    
    # Generate flow channels
    for i in range(start_pos+ spacing//2, size, spacing):
        if i + width > size:
            break
        grid[:, i:i+width] = 0
    
    # Enforce boundaries as obstacles
    #grid[0, :] = 1; grid[-1, :] = 1; grid[:, 0] = 1; grid[:, -1] = 1
    return grid

straight_grid = generate_straight_channels(size, width, spacing)
#straight_grid = np.rot90(straight_grid)

# Add inlet and output transition layers
zeros_row = np.zeros((5,size))
straight_grid = np.vstack([zeros_row,straight_grid,zeros_row])
straight_grid = np.rot90(straight_grid)

plt.figure(figsize=(8, 8))
straight_grid = np.ma.masked_where((straight_grid > 0.0), straight_grid)
plt.imshow(straight_grid, vmin= -1.0,
                    vmax= 1.0, aspect='auto', interpolation='nearest')
#plt.show()

filename = 'straight_'+str(size)+'_'+str(spacing)
plt.axis('off')
plt.savefig(filename+'_rot90'+'.png', dpi=300 ,bbox_inches='tight', pad_inches=0)
#plt.show()

# Flattening matrix (row priority)
flat_straight = straight_grid.flatten()

# Save as txt file
np.savetxt(filename+'_rot90'+'.txt', flat_straight, fmt="%d") 
