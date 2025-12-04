# Single-phase-LBM-for-gas-seepage-in-porous-media
This is a single-phase LBM code suitable for gas seepage in porous media

Reference code sources:
- https://github.com/jviquerat/lbm/workflows/lbm/badge.svg?branch=master
- https://github.com/Jeff-Hugh/GenPorMed.git

Add several functions based on the reference codes

## Contents
Cases implementations include:
- Couette flow (../app/huayi.py)
- Poiseuille flow (../app/poiseuille.py)
- Around circular cylinder flow (../app/turek.py)
- Variable cross-section slit flow (../app/step.py)

Something newly added:

- View LBM principle (unit.py)
- Influence of Knudsen number and gas property parameters (../app/knudsen.py)
- Complex porous media without slip effect (../app/porous.py & porous_1.py)
- Complex porous media with slip effect (../app/porous_kn.py)

## Code architecture
- Program entry (start.py & ../core/run.py)
- Case setting (../app)
- Grid import (../app/base_app.py)
- Solver (../core)
- Matrix algorithm (../core)
- Output (../plot/plot.py)
- Simple shape drawing (../utils)
- Complex pore structure generation

## How to run LBM simulation
1. Clone the repo
```shell
    git clone https://github.com/whupmcww/Single-phase-LBM-for-gas-seepage-in-porous-media.git
```

2. Select which case you want to use. Take case porous.py as example, set parameters for your case. Specify the grid file path if available
```shell
    class porous(base_app):
    def __init__(self):
        ...
        # Set parameters for YourCase
        self.name     = 'porous'
        self.L_real   = 1.0e-2
        self.L_lbm    = 1000

        # Obstacles
        self.filepath = 'path to your grid file'
        self.filename = self.filepath +'YourGridfile.txt'
        self.obstacle_ww =  obstacle('users_define')
        ...
```

3. Modify the case name in the start.py
```shell
    if __name__ == '__main__':
        ...
        # Instanciate app
        app = app_factory.create('porous')
        ...
```

4. Run the file by IDLE F5 or Powershell command. View the result file output in the specified folder
```shell
    python start.py
```

## How to generate porous grid by QSGS
1. Use code in test_gen
2. Set the generation mode and QSGS parameters in run.cpp
```shell
    ...
    void Generate2D()
    {
        const int M = 50;
	    const int N = 50;
	    const double phi = 0.371;			/// target porosity
	    const double p_cd = 0.05;		/// core distribution probability
        ...
        Porous2D porous(M, N, phi, p_cd, z_h, z_v, f);
        ...
        char filename[100] = "./output/Grid.txt";
	    porous.output2datatxt(M, N, s, filename);
    }
    ...
```
3. Compile code using Cmake and run, you can put these command in a .bat file
```shell
    cmake ./
    make
    ./run
```

## Case explanation
Complex porous pressure and velocity cloud map
<p align="center">
    <img width="300" alt="Pressure" src="porous_kn_results/Directionality/p_norm_622.png"> <img width="300" alt="Velocity" src="porous_kn_results/Directionality/u_norm_622.png"> <img width="300" alt="Permeability" src="porous_kn_results/Directionality/k_622000.png">
</p>
