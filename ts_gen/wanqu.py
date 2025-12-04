import numpy as np
import matplotlib.pyplot as plt

size=1000
width=10
spacing=30
amplitude=10
frequency=0.01

def generate_curved_channels(size, width, spacing, amplitude, frequency):
    """生成正弦弯曲流道，与直流道保持相同数量"""
    grid = np.ones((size, size))
    
    # 使用与直流道相同的布局参数
    channel_count = (size - width) // spacing
    start_pos = (size - (channel_count * spacing + width)) // 2
    
    for row in range(size):
        # 正弦波参数（一个完整周期）
        sin_offset = int(amplitude * np.sin(row * 2*np.pi/(size/16)))
        
        # 生成交错流道（位置在直流道之间）
        for i in range(start_pos + spacing//2, size, spacing):
            left = max(0, i + sin_offset - width//2)
            right = min(size, i + sin_offset + width//2)
            
            # 确保最小有效宽度
            if right - left >= width//2:
                grid[row, left:right] = 0
    
    # 边界保护
    #grid[0, :] = 1; grid[-1, :] = 1; grid[:, 0] = 1; grid[:, -1] = 1
    return grid

curved_grid = generate_curved_channels(size, width, spacing, amplitude, frequency)
#curved_grid = np.rot90(curved_grid)

#添加边界层
zeros_row = np.zeros((5,size))
curved_grid = np.vstack([zeros_row,curved_grid,zeros_row])
curved_grid = np.rot90(curved_grid)

plt.figure(figsize=(8, 8))
curved_grid = np.ma.masked_where((curved_grid > 0.0), curved_grid)
plt.imshow(curved_grid, vmin= -1.0,
                    vmax= 1.0, aspect='auto', interpolation='nearest')


filename = 'curved_'+str(size)+'_'+str(spacing)+'_16'
plt.axis('off')
plt.savefig(filename+'_rot90'+'.png', dpi=300 ,bbox_inches='tight', pad_inches=0)
#plt.show()

# 展平矩阵（按行优先）
flat_curved = curved_grid.flatten()

# 保存为单列 txt 文件
np.savetxt(filename+'_rot90'+'.txt', flat_curved, fmt="%d")  # %d 表示整数格式

