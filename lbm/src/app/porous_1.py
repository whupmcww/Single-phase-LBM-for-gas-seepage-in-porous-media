# Generic imports
import math

# Custom imports
from lbm.src.app.base_app  import *
from lbm.src.core.lattice_obstacle  import *
from lbm.src.core.obstacle_ww import *
from lbm.src.utils.buff    import *
from lbm.src.plot.plot     import *

###############################################
### A step with modified inlets/outlets
class porous_1(base_app):
    def __init__(self):

        # Free arguments
        self.name        = 'porous_1'
        self.L_real      = 1.0e-2
        self.L_lbm       = 1000
        self.Re_lbm      = 0.1
        self.t_max       = 1.0 # Real cacu time
        self.nu_real     = 1.06e-4 # Real fluid viscosity
        self.ny          = self.L_lbm 
        self.dy          = 1.0
        self.rho_lbm     = 1.0 # LBM density
        self.rho_real    = 1.0e3 # Real density
        self.Cs_real     = 1483
        self.tau_lbm     = 0.93
        
        self.stop        = 'it'
        self.obs_cv_ct   = 1.0e-3
        self.obs_cv_nb   = 1000

        # Output parameters
        self.output_freq = 1000
        self.output_it   = 0
        self.dpi         = 200
        self.cacu_step = []
        self.output_var = []

        # Deduce remaining lbm parameters
        self.compute_lbm_parameters()

        # Obstacles
        self.filepath = 'G:/LBM/para_test/sample/k_test/kuozhan/'
        self.filename = self.filepath +'100_100_0.5_0.01_0.1_0.01_0.025_10x10'+'.txt'
        self.obstacle_ww =  obstacle('users_define')

    ### Compute remaining lbm parameters
    def compute_lbm_parameters(self):
        self.Cs      = 1.0/math.sqrt(3.0)
        self.Cx      = self.L_real/self.L_lbm
        self.nu_lbm  = self.Cs**2*(self.tau_lbm-0.5)
        self.Cnu     = self.nu_real/self.nu_lbm
        self.Ct      = self.Cx**2/self.Cnu
        self.dt      = self.Ct
        #self.u_lbm   = self.Re_lbm*self.nu_lbm/self.L_lbm
        self.Cu      = self.Cx/self.Ct
        #self.u_real  = self.u_lbm*self.Cu
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
        lattice.u[:,np.where(lattice.lattice > 0.0)] = 0.0
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

        ############## Velocity inlet
        #lattice.u_left[0,:] = ret*self.u_lbm
        #lattice.u_left[1,:] = 0.0
        ############## Pressure inlet
        self.p_left = 10.0e6
        lattice.rho_left[:] = self.p_left / (self.rho_real*self.Cs_real**2) /1000.0 + self.rho_lbm
        #lattice.p_left[:] = lattice.rho_left[:] * self.Cs**2
        #print(lattice.rho_left[1])
        lattice.u_left[1,:] = 0.0

        lattice.u_top[0,:]   = 0.0
        lattice.u_bot[0,:]   = 0.0
        lattice.u_right[1,:] = 0.0
        lattice.rho_right[:] = self.rho_lbm

    ### Set boundary conditions
    def set_bc(self, lattice):

        # Obstacle
        lattice.bounce_back_obstacle_ww(self.obstacle_ww)

        # Wall BCs
        lattice.zou_he_bottom_wall_velocity()

        ############### Velocity inlet
        #lattice.zou_he_left_wall_velocity()
        ############### Pressure inlet
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
        export_poro_flux(lattice, it, self.cacu_step, self.output_var, self.nu_lbm, self.output_freq)

        # Increment plotting counter
        self.output_it += 1

