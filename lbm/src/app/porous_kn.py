# Generic imports
import math

# Custom imports
from lbm.src.app.base_app  import *
from lbm.src.core.lattice_obstacle_addForce  import *
from lbm.src.core.obstacle_ww import *
from lbm.src.utils.buff    import *
from lbm.src.plot.plot     import *

###############################################
### A step with modified inlets/outlets
class porous_kn(base_app):
    def __init__(self):

        # Free arguments
        self.name        = 'porous_kn'
        self.L_real      = 0.1e-5
        self.L_lbm       = 1000
        self.Re_lbm      = 0.1
        self.t_max       = 1.0 #实际计算时间
        self.nu_real     = 1.06e-4 #实际流体黏度
        self.ny          = self.L_lbm 
        self.dy          = 1.0
        self.rho_lbm     = 1.0 #LBM密度
        self.rho_real    = 0.08988 #1.0e3 #实际密度
        self.Cs_real     = 340
        #self.tau_lbm     = 0.93

        ########### 计算松弛时间所需要的参数(仅对气体)
        self.m_gas       = 3.35E-27
        self.d_gas       = 0.29E-9
        self.N           = 10

        self.T = 300.0
        self.R = 8.314
        self.Tc = 33.18
        self.pc = 1.315e6
        self.omega = -0.128
        ###########
        
        #self.stop        = 'it'
        self.stop        = 'k'
        self.k_setError  = 1.0e-4
        self.k_error     = 1.0
        self.k_error_1   = []
        self.it_error_1  = []
        self.judge       = True
        self.obs_cv_ct   = 1.0e-3
        self.obs_cv_nb   = 1000

        # Output parameters
        self.output_freq = 5000
        self.output_it   = 0
        self.dpi         = 200
        self.cacu_step = []
        self.output = []

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

        # Obstacles
        self.filepath = 'G:/LBM/para_test/sample/渗透率测试/kuozhan/方向性/'
        self.filename = self.filepath +'100_100_0.3_0.01_0.1_0.001_0.025_10x10_new'+'.txt'
        self.obstacle_ww =  obstacle('users_define')

    ### Compute remaining lbm parameters
    def compute_lbm_parameters(self):
        ############ 计算松弛时间(仅对气体)
        self.H           = self.L_real / 1000 * 10
        self.b_gas       = 2.0*math.pi*self.d_gas**3/(3.0*self.m_gas)
        self.b_rho       = self.b_gas*self.rho_real
        self.modify      = 1.0+5.0/8.0*self.b_rho+0.2869*self.b_rho**2+0.1103*self.b_rho**3+0.0386*self.b_rho**4
        self.lamda       = self.m_gas/(math.sqrt(2)*self.rho_real*self.d_gas**2*self.modify)
        self.k_n         = self.lamda/self.H
        print(self.k_n)
        self.fi_kn       = 1.0 / (1.0 + 2.0 * self.k_n)
        self.tau_lbm     = 0.5 + math.sqrt(6.0/math.pi)*self.k_n*self.fi_kn*self.N*(1.0+0.5*self.b_rho*self.modify)**2/self.modify
        #print(self.tau_lbm)
        #self.tau_lbm     = 0.93

        self.const   = 0.0778*self.R*self.Tc/self.pc
        self.alpha  = (1+(0.37464+1.54226*self.omega-0.26992*self.omega**2)+(1.0-(self.T/self.Tc)**0.5))**2
        self.a = 0.45724*self.R**2*self.Tc**2/self.pc
        ############

        ############ 计算r
        self.first = (1.0/self.N)**2/(4.0*self.k_n*self.fi_kn)
        self.sigma = 1.0
        self.sigma_v = (2.0-self.sigma)/self.sigma
        self.B1 = (1.0-0.1817*self.sigma)
        self.second = self.B1*self.sigma_v
        self.B2 = 0.8
        self.third = (2.0*self.B2-8.0*(1.0+0.5*self.b_rho*self.modify)**4/(math.pi*self.modify**2))*self.k_n*self.fi_kn
        self.r_0 = 1.0+math.sqrt(math.pi/6.0)*self.modify/(1.0+0.5*self.b_rho*self.modify**2*(self.first+self.second+self.third))
        self.r = 1.0/self.r_0
        #self.r = 1.0
        ############

        self.Cs      = 1.0/math.sqrt(3.0)
        self.Cx      = self.L_real/self.L_lbm
        self.nu_lbm  = self.Cs**2*(self.tau_lbm-0.5)
        self.Cnu     = self.nu_real/self.nu_lbm
        self.Ct      = self.Cx**2/self.Cnu
        self.dt      = self.Ct
        #self.u_lbm   = self.Re_lbm*self.nu_lbm/self.L_lbm
        self.Cu      = self.Cx/self.Ct
        #self.u_real  = self.u_lbm*self.Cu
        self.Crho    = self.rho_real / self.rho_lbm
        self.nx      = self.ny
        self.dx      = self.dy
        self.it_max  = math.floor(self.t_max/self.dt)
        self.sigma   = math.floor(10*self.nx)

    ### Add obstacles and initialize fields
    def initialize(self, lattice):

        # Add obstacles to lattice
        self.add_obstacles_load(lattice, self.filename,self.obstacle_ww)

        # Initialize fields
        #self.set_inlets(lattice, 0)
        mask = lattice.lattice > 0.0
        lattice.u[:, mask] = 0.0
        #lattice.u[:,np.where(lattice.lattice > 0.0)] = 0.0
        lattice.rho *= self.rho_lbm
        self.set_inlets(lattice, 0)

        # Output image
        lattice.generate_image(self.obstacle_ww)

        # Compute first equilibrium
        lattice.equilibrium()
        lattice.g = lattice.g_eq.copy()

    ### Set inlet fields
    def set_inlets(self, lattice, it):

        lx = lattice.lx
        ly = lattice.ly

        val  = it
        ret  = (1.0 - math.exp(-val**2/(2.0*self.sigma**2)))

        ############## 速度入口
        #lattice.u_left[0,:] = ret*self.u_lbm
        #lattice.u_left[1,:] = 0.0
        ##############

        ############## 压力入口
        self.p_left = 1.0e6
        lattice.rho_left[:] = self.p_left / (self.rho_real*self.Cs_real**2) /1000.0 + self.rho_lbm
        #lattice.p_left[:] = lattice.rho_left[:] * self.Cs**2
        #print(lattice.rho_left[1])
        lattice.u_left[1,:] = 0.0
        ##############

        lattice.u_top[0,:]   = 0.0
        lattice.u_bot[0,:]   = 0.0
        lattice.u_right[1,:] = 0.0
        lattice.rho_right[:] = self.rho_lbm

    ### Set boundary conditions
    def set_bc(self, lattice):

        # Obstacle
        #lattice.bounce_back_obstacle_ww(self.obstacle_ww)
        lattice.bounce_back_obstacle_kn_ww(self.obstacle_ww)

        # Wall BCs
        lattice.zou_he_bottom_wall_velocity()

        ###############速度入口
        #lattice.zou_he_left_wall_velocity()
        ###############压力入口
        lattice.zou_he_left_wall_pressure()

        lattice.zou_he_top_wall_velocity()
        lattice.zou_he_right_wall_pressure()
        lattice.zou_he_bottom_left_corner()
        lattice.zou_he_top_left_corner()
        lattice.zou_he_top_right_corner()
        lattice.zou_he_bottom_right_corner()

    ### Write outputs
    def outputs(self, lattice, it):

        # Check iteration
        if (it%self.output_freq != 0): return

        # Output field
        plot_norm_p(lattice, 0.0, 1.5, self.output_it, self.dpi)
        plot_norm_u(lattice, 0.0, 1.5, self.output_it, self.dpi)

        output_it, output_var = lattice.cacu_k(it)
        self.cacu_step.append(output_it)
        self.output.append(output_var)

        export_poro_flux(lattice, it, self.cacu_step, self.output, self.output_freq)

        # Increment plotting counter
        self.output_it += 1

    ### 误差判断
    def judge_kError(self, it):
        k_error = self.cacu_kError()
        if k_error <= self.k_setError:
            self.k_error_1.append(k_error)
            self.it_error_1.append(it)
            if it == self.it_error_1[0]+10000:
                if self.k_error_1[len(self.k_error_1)-1] <= self.k_setError:
                    self.judge = False

        return self.judge
    
    ### Caculate k_error
    def cacu_kError(self):
        n = len(self.cacu_step)
        if n >= 2:   
            last_k = self.output[n-2]
            now_k =  self.output[n-1] 
            if now_k != 0:  
                self.k_error = math.fabs((now_k - last_k)/ now_k)

        return self.k_error
    
    def export_k(self, lattice):
        export_k_as_csv(lattice, self.cacu_step,self.output)


