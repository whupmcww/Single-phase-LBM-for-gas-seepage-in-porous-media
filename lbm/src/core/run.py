# Generic imports
import os
import time

# Custom imports
from lbm.src.app.app      import *
# from lbm.src.core.lattice_addForce import *
# from lbm.src.core.lattice_obstacle_addForce import *

########################
# Run lbm simulation
########################
def run(lattice, app):

    # Initialize fields and distributions
    app.initialize(lattice)

    # Timer and loop data
    start_time = time.time()
    it         = 0
    compute    = True

    # Solve
    print('### Solving')
    try:    
        while (compute):
            # Printings
            #app.printings(it)

            # Set inlets
            app.set_inlets(lattice, it)

            # Compute force
            #lattice.cacu_force()
            #print(lattice.Force)

            # Compute macroscopic fields
            lattice.macro()
            #print(lattice.u[0])

            # Compute permeability
            lattice.cacu_k(it)

            # Output field
            app.outputs(lattice, it)

            # Judge error for permeability k
            app.judge_kError(it)

            # Compute equilibrium state
            lattice.equilibrium()

            # Streaming
            lattice.collision_stream()

            # Boundary conditions
            app.set_bc(lattice)

            # Compute observables (drag, lift, etc)
            #app.observables(lattice, it)

            # Check stopping criterion
            compute = app.check_stop(it)


            if not compute:
                app.export_k(lattice)

            # Increment iteration
            it += 1

    except KeyboardInterrupt:
        app.export_k(lattice)

    finally:
        # Count time
        end_time = time.time()
        print("# Loop time = {:f}".format(end_time - start_time))

        # Perform final operations and outputs
        app.finalize(lattice)



