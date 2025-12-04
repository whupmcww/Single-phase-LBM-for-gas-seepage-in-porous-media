# Generic imports
import os
import math
import csv
import numpy              as np
import matplotlib.pyplot  as plt

kuozhan = 10

filepath='G:/LBM/para_test/sample/渗透率测试/方向性/50x50/rand=15/' #输入路径
filepath1='G:/LBM/para_test/sample/渗透率测试/' #输出路径
qianzhui= '50_50_0.3_0.05_0.1_0.01_0.025'
name = filepath+qianzhui
filename = name + '.txt'
name_1 = filepath1+'kuozhan/方向性/50x50/rand=15/'+ qianzhui
name_2 = name_1+'_'+str(kuozhan)+'x'+str(kuozhan)#+'_rot90'
name_3 = name_2+'.txt'

### ************************************************
### Output 2D flow amplitude
def plot_mesh():
    data_size_x = 50
    data_size_y = 50

    in_dat = np.loadtxt(filename)
    in_dat[in_dat>0] = 1
    in_dat = np.reshape(in_dat, (data_size_x,data_size_y),order='F')
    in_dat = np.rot90(in_dat)

    #扩展矩阵
    in_dat = np.kron(in_dat, np.ones((kuozhan, kuozhan)))

    #添加边界层
    zeros_row = np.zeros((5, data_size_x*kuozhan))
    in_dat = np.vstack([zeros_row, in_dat, zeros_row])
    in_dat = np.rot90(in_dat)

    #重新写出
    flattened_data = in_dat.flatten(order='C')
    with open(name_3, 'w') as file:
        for num in flattened_data:
            file.write(f"{num}\n")
   
    # Mask obstacles
    data = np.ma.masked_where((in_dat > 0.0), in_dat)
    data = np.rot90(data)
    #print(data)

    # Plot
    plt.clf()
    
    image_height, image_width = data.shape
    fig, ax = plt.subplots(figsize=(image_width/10.0, image_height/10.0))

    ax.imshow(data,cmap='RdBu_r', aspect='auto')

    
    #ax.grid(True, which='both', linestyle='-', linewidth=0.5)
    ax.grid(False, which='both')
    ax.tick_params(which='both', top=False, right=False, bottom=False, left=False,
            labelbottom=False, labelleft=False)
    ax.tick_params(which='both', top=False, right=False, bottom=False, left=False,
            grid_color='black', grid_linestyle='-', grid_linewidth=0.5,labelbottom=False, labelleft=False)
    ax.set_xticks(np.arange(0.5, data_size_x*kuozhan+2, 1.01))
    ax.set_yticks(np.arange(0.5, data_size_y*kuozhan, 1.0))

    #plt.axis('off')
    
    filename1 = name_2+'.png'
    #plt.savefig(filename1, dpi=300)
    plt.savefig(filename1)
    plt.pause(0.5)
    plt.close()
    print('生成完成')

plot_mesh()

def gaixie():
    # 打开文件并读取内容
    with open(name_3, 'r', encoding='utf-8') as file:
        content = file.read()

    # 替换内容
    content = content.replace('0.0', '0').replace('1.0', '1')

    # 将修改后的内容写回文件
    with open(name_2+'_new'+'.txt', 'w', encoding='utf-8') as file:
        file.write(content)

    print('改写完成')

gaixie()



