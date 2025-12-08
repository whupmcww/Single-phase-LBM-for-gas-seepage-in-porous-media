# Generic imports
import os
import math
import csv
import numpy              as np
import matplotlib.pyplot  as plt
import matplotlib.ticker as ticker

### ************************************************
### Output 2D flow amplitude
def plot_norm_p(lattice, val_min, val_max, output_it, dpi):

    # Compute norm
    v = np.sqrt(lattice.u[0,:,:]**2+lattice.u[1,:,:]**2)
    p = (lattice.rho[:,:]-lattice.rho_lbm)*lattice.Cs_real**2 * 1000.0 * lattice.rho_real / 1e6

    # Mask obstacles
    #v[np.where(lattice.lattice > 0.0)] = -1.0
    vm = np.ma.masked_where((lattice.lattice > 0.0), v)
    vm = np.rot90(vm)
    
    p = np.ma.masked_where((lattice.lattice > 0.0), p)
    p = np.rot90(p)

    # Plot
    plt.clf()
    # fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    # fig.subplots_adjust(0,0,1,1)
    image_height, image_width = vm.shape
    fig, ax = plt.subplots(figsize=(image_width/ 40, image_height/ 50))
    #fig, ax = plt.subplots(figsize=(image_width, image_height))
    im=ax.imshow(p,cmap='RdBu_r', aspect='auto')
            #     cmap='RdBu_r',
            #    interpolation = 'spline16')
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label('P Values / MPa', fontsize=24)
    # def format_func(value, tick_number):
    #     return f'{value:.2e}'
    # cbar.formatter = ticker.FuncFormatter(format_func)
    
    #ax.grid(True, which='both', linestyle='-', linewidth=0.5)
    ax.grid(False, which='both')
    # ax.tick_params(which='both', top=False, right=False, bottom=False, left=False,
    #         labelbottom=False, labelleft=False)
    # ax.tick_params(which='both', top=False, right=False, bottom=False, left=False,
    #         grid_color='black', grid_linestyle='-', grid_linewidth=0.5,labelbottom=False, labelleft=False)
    # ax.set_xticks(np.arange(0.5, lattice.nx, 1.0))
    # ax.set_yticks(np.arange(0.5, lattice.ny, 1.0))

    filename = lattice.png_dir+'p_norm_'+str(output_it)+'.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi ,bbox_inches='tight', pad_inches=0)
    #plt.pause(0.5)
    plt.close()

### ************************************************
### Output 2D flow amplitude
def plot_norm_u(lattice, val_min, val_max, output_it, dpi):

    # Compute norm
    v = np.sqrt(lattice.u[0,:,:]**2+lattice.u[1,:,:]**2)

    # Mask obstacles
    #v[np.where(lattice.lattice > 0.0)] = -1.0
    vm = np.ma.masked_where((lattice.lattice > 0.0), v)
    vm = np.rot90(vm)

    # Plot
    plt.clf()
    # fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    # fig.subplots_adjust(0,0,1,1)
    image_height, image_width = vm.shape
    fig, ax = plt.subplots(figsize=(image_width/ 40, image_height/ 50))
    #fig, ax = plt.subplots(figsize=(image_width, image_height))
    im=ax.imshow(vm,cmap='RdBu_r', aspect='auto')
            #     cmap='RdBu_r',
            #    interpolation = 'spline16')
    cbar = plt.colorbar(im, ax=ax)
    cbar.ax.tick_params(labelsize=18)
    cbar.set_label('u Values / m/s', fontsize=24)
    # def format_func(value, tick_number):
    #     return f'{value:.2e}'
    # cbar.formatter = ticker.FuncFormatter(format_func)
    
    #ax.grid(True, which='both', linestyle='-', linewidth=0.5)
    ax.grid(False, which='both')
    # ax.tick_params(which='both', top=False, right=False, bottom=False, left=False,
    #         labelbottom=False, labelleft=False)
    # ax.tick_params(which='both', top=False, right=False, bottom=False, left=False,
    #         grid_color='black', grid_linestyle='-', grid_linewidth=0.5,labelbottom=False, labelleft=False)
    # ax.set_xticks(np.arange(0.5, lattice.nx, 1.0))
    # ax.set_yticks(np.arange(0.5, lattice.ny, 1.0))

    filename = lattice.png_dir+'u_norm_'+str(output_it)+'.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi ,bbox_inches='tight', pad_inches=0)
    #plt.pause(0.5)
    plt.close()

### ************************************************
### Output 2D flow contour
def plot_contour(lattice, output_it, dpi):
    
    x  = np.linspace(0, 1, lattice.nx)
    y  = np.linspace(0, 1, lattice.ny)
    ux = lattice.u[0,:,:].copy()
    uy = lattice.u[1,:,:].copy()
    uy = np.rot90(uy)
    ux = np.rot90(ux)
    uy = np.flipud(uy)
    ux = np.flipud(ux)
    vm = np.sqrt(ux**2+uy**2)
    rho = lattice.rho[:,:].copy()
    rho = np.rot90(rho)
    rho = np.flipud(rho)

    # Plot
    plt.clf()
    fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    fig.subplots_adjust(0,0,1,1)

    plt.contour(x, y, rho, cmap='RdBu_r',
               )
    filename = lattice.png_dir+'u_ctr_'+str(output_it)+'.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi)
    plt.close()

### ************************************************
### Output 2D streamlines
def plot_streamlines(lattice, dpi):

    # The outputted streamplot is rotated and flipped...
    plt.clf()
    fig, ax = plt.subplots(figsize=plt.figaspect(vm))
    fig.subplots_adjust(0,0,1,1)
    ux = lattice.u[0,:,:].copy()
    uy = lattice.u[1,:,:].copy()
    uy = np.rot90(uy)
    ux = np.rot90(ux)
    uy = np.flipud(uy)
    ux = np.flipud(ux)
    vm = np.sqrt(ux**2+uy**2)
    vm = np.rot90(vm)
    x  = np.linspace(0, 1, lattice.nx)
    y  = np.linspace(0, 1, lattice.ny)
    nn = 20
    u  = np.linspace(0, 1, nn)
    str_pts = []
    for a in range(nn):
        for b in range(nn):
            str_pts.append([u[a],u[b]])
            plt.streamplot(x, y, ux, uy,
                           linewidth    = 0.2,
                           cmap         = 'RdBu_r',
                           arrowstyle   = '-',
                           start_points = str_pts,
                           density      = 20)

    filename = lattice.output_dir+'u_stream.png'
    plt.axis('off')
    plt.savefig(filename, dpi=dpi)
    plt.close()


def plot_outlet_velocity_curve(lattice, dpi, it):
    half = math.floor(0.5*lattice.lx)
    ux = lattice.u[0,half, :].copy()
    y  = np.linspace(0, 1, lattice.ny)

    # p_left = lattice.rho[0, :].sum() / lattice.ny
    # p_right = lattice.rho[lattice.lx, :].sum() / lattice.ny
    # # det_p = ratio_p2rho*(p_left - p_right) * (lattice.Cs**2)
    # # print(p_left - p_right)
    # valid_ux = []
    # v_ux = 0
    # ratio_p2rho = 1.0
    # for i in range(lattice.ny):
    #     det_p = ratio_p2rho*(p_left - p_right) * (lattice.Cs**2)
    #     v_ux = det_p*(y[i]*lattice.ny)*(lattice.ny-(y[i]*lattice.ny))/(2.0*lattice.nu_lbm)
    #     #valid_ux.append(v_ux)
    #     if lattice.ny%2==0:
    #         if i == 0.5*lattice.ny-1:
    #             if ux[i]:
    #                 ratio_p2rho = v_ux / ux[i]
    #     else:
    #         if i == 0.5*(lattice.ny-1):
    #             if ux[i]:
    #                 ratio_p2rho = v_ux / ux[i]
    # #print(ratio_p2rho)
    # for i in range(lattice.ny):
    #     det_p = (p_left - p_right) * (lattice.Cs**2)/ratio_p2rho
    #     v_ux = det_p*(y[i]*lattice.ny)*(lattice.ny-(y[i]*lattice.ny))/(2.0*lattice.nu_lbm)
    #     valid_ux.append(v_ux)
    
    # Plot
    plt.clf()
    fig, ax = plt.subplots()
    #fig.subplots_adjust(0,0,1,1)
    #color_cycle = plt.rcParams['axes.prop_cycle']
    #for color in color_cycle:
        #print(color['color'])

    # ax.plot(y, valid_ux, linestyle='', marker='x', markersize=8, markeredgewidth=2, color='#ff7f0e', label='Theory')
    ax.plot(y, ux, color='#1f77b4', label='LBM')
    ax.legend()
    ax.set_xlabel('ux')
    ax.set_ylabel('y/H')
    
    #plt.show()
    filename = lattice.png_dir+'kn'+'_'+str(lattice.k_n)+'_'+str(it)+'.png'
    plt.axis('on')
    plt.savefig(filename, dpi=dpi)
    plt.close()

    filename1 = lattice.png_dir+'kn'+'_'+str(lattice.k_n)+'_'+str(it)+'.csv'
    with open(filename1, mode='w', newline='') as file:
        writer = csv.writer(file)
    
        # writer.writerow(['x', 'y'])
        # for x1, y1, y2 in zip(y, ux, valid_ux):
        #     writer.writerow([x1, y1, y2])
        for x1, y1 in zip(ux, y):
              writer.writerow([x1, y1])

    print("CSV file successfully written")

################################################
def plot_outlet_pressure_curve(lattice, dpi, it):
    half = math.floor(0.5*lattice.ly)
    rho = lattice.rho[:, half].copy()
    x  = np.linspace(0, 1, lattice.nx)
    
    # Plot
    plt.clf()
    fig, ax = plt.subplots()
    #fig.subplots_adjust(0,0,1,1)
    #color_cycle = plt.rcParams['axes.prop_cycle']
    #for color in color_cycle:
        #print(color['color'])

    ax.plot(x, rho, color='#1f77b4', label='LBM')
    ax.legend()
    ax.set_xlabel('p')
    ax.set_ylabel('x/L')
    
    #plt.show()
    filename = lattice.png_dir+'p'+'_'+str(it)+'_'+'.png'
    plt.axis('on')
    plt.savefig(filename, dpi=dpi)
    plt.close()

    filename1 = lattice.png_dir+'p'+str(lattice.u_lbm)+'_'+str(it)+'_'+'.csv'
    with open(filename1, mode='w', newline='') as file:
        writer = csv.writer(file)
        for x1, y1 in zip(x, rho):
            writer.writerow([x1, y1])

    print("CSV file successfully written")

################################################
def plot_cavity_validation(lattice, dpi, it):
    half = math.floor(0.5*lattice.lx)
    ux = lattice.u[0,half, :].copy()
    ux_1 = ux/ux[lattice.ly]
    y  = np.linspace(0, 1, lattice.ny)
    #validation_value
    valid_ux = []
    v_ux = 0
    for i in range(lattice.ny):
        b = lattice.ly / 4.5
        k0 = 1.0/(1+lattice.ly/b)
        #v_ux = (y[i]*ux[lattice.ly]*(1-k0)+k0*ux[lattice.ly])/ux[lattice.ly]
        v_ux = y[i]
        valid_ux.append(v_ux)
        
    
    # Plot
    plt.clf()
    fig, ax = plt.subplots()
    #fig.subplots_adjust(0,0,1,1)
    #color_cycle = plt.rcParams['axes.prop_cycle']
    #for color in color_cycle:
        #print(color['color'])

    ax.plot(valid_ux, y, linestyle='', marker='x', markersize=8, markeredgewidth=2, color='#ff7f0e', label='Theory')
    ax.plot(ux_1, y,color='#1f77b4', label='LBM')
    ax.legend()
    ax.set_xlabel('ux')
    ax.set_ylabel('y/H')
    
    #plt.show()
    filename = lattice.png_dir+'u_out'+str(lattice.u_lbm)+'_'+str(it)+'_'+'.png'
    plt.axis('on')
    plt.savefig(filename, dpi=dpi)
    plt.close()

    filename1 = lattice.png_dir+'u_out'+str(lattice.u_lbm)+'_'+str(it)+'_'+'.csv'
    with open(filename1, mode='w', newline='') as file:
        writer = csv.writer(file)

        # writer.writerow(['x', 'y'])
        for x1, y1, y2 in zip(y, ux_1, valid_ux):
            writer.writerow([x1, y1, y2])

    print("CSV file successfully written")

################################################
def export_g(lattice, it):
    for i in range(9):
        # transposed_array = list(zip(*lattice.g[i]))
        # transposed_array = [list(row) for row in transposed_array]

        filename1 = lattice.png_dir+'g_out'+'['+str(i)+']'+'_'+str(it)+'_'+'.csv'
        with open(filename1, mode='w', newline='') as file:
            writer = csv.writer(file)
            
            # for row in transposed_array:
            #     writer.writerow(row)
            for row in lattice.g[i]:
                writer.writerow(row)

    print("CSV file successfully writte")

################################################
def export_step_velocity(lattice, it):
    if it%200==0:
        half_1 = math.floor(1.0/16.0*lattice.lx)
        half_2 = math.floor(2.5/16.0*lattice.lx)
        half_3 = math.floor(5.0/16.0*lattice.lx)
        ux_1 = lattice.u[0,half_1, :].copy()
        ux_2 = lattice.u[0,half_2, :].copy()
        ux_3 = lattice.u[0,half_3, :].copy()
        y =  np.linspace(0, 1, lattice.ny)

        assert len(y) == len(ux_1) == len(ux_2) == len(ux_3)

        data = [
            ["y", "ux_1", "ux_2", "ux_3"],  # Lables
            *[list(row) for row in zip(y, ux_1, ux_2, ux_3)]  # Data rows
            ]

        filename1 = lattice.output_dir+'v_out'+str(it)+'.csv'
        with open(filename1, mode='w', newline='') as file:
            writer = csv.writer(file)
                
            writer.writerows(data)

        print("CSV file successfully writte")
    
###################################################
def export_poro_flux(lattice, it, cacu_step, output, plot_every):
    # out = math.floor(lattice.lx)
    # ux = lattice.u[0,out, :].copy() # Inlet velocity
    # rho_inlet = lattice.rho[0,:].copy()
    # rho_outlet = lattice.rho[out,:].copy()
    # y  = np.linspace(0, 1, lattice.ny)

    # out_area = 0
    # u_total = 0.0
    # rho_in_total = 0.0
    # rho_out_total = 0.0
    # for j in range(lattice.ly):
    #     if ux[j] != 0:
    #         out_area += 1
    #         u_total += ux[j]
    #         rho_in_total += rho_inlet[j]
    #         rho_out_total += rho_outlet[j]
    
    # u_avg = u_total / lattice.ly
    # rho_in_avg = rho_in_total
    # rho_out_avg = rho_out_total
    # out_flux = u_total
    # det_rho = rho_in_avg - rho_out_avg
    
    # if det_rho ==0:
    #     k = 0
    # else:
    #     k = out_flux / lattice.ly * lattice.lx * nu / (3.0*det_rho)

    # output_var = k

    # cacu_step.append(it)
    # output.append(output_var)
    if(it % plot_every == 0):
        fig, ax1 = plt.subplots()
        ax1.plot(cacu_step,output , marker='o', linestyle='-', color='b', label='A')
        filename = lattice.png_dir+'k_'+str(it)+'.png'
        plt.savefig(filename, dpi=200, bbox_inches='tight', pad_inches=0)
        #plt.pause(0.5)
        plt.close()

        
###########################################
def export_k_as_csv(lattice, step, output):

    data = [
            ["step", "k"],  # Lables
            *[list(row) for row in zip(step, output)]  # Data rows
            ]

    filename = lattice.output_dir+'k_out'+'.csv'
    with open(filename, mode='w', newline='') as file:
        writer = csv.writer(file)
                
        writer.writerows(data)

    print("CSV file successfully writte")
