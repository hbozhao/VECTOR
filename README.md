# VECTOR
Phase field simulation codes for "Condensate interfacial forces reposition DNA loci and probe chromatin viscoelasticity" by Amy R. Strom*, Yoonji Kim*, Hongbo Zhao, Yi-Che Chang, Natalia D. Orlovsky, Andrej Košmrlj, Cornelis Storm, Clifford P. Brangwynne

* To perform simulations using different viscoelasticity models, use the script `VECTOR_study.m`.
  In `VECTOR_study.m`, users may run the base case using the Newtonian model by setting `run_base_case = true`. The base case simulation consists of three steps:
  1. equilibration from the initial condition when the light is off
  2. equilibriation with the light on, until the Corelet droplet forms and merges into one
  3. turn off the light.

* To run other other viscoelastic models, set `run_base_case = false` and set `study` to the name of the viscoelastic model (see `VECTOR_study.m` for details). These studies use the end of the second step as the initial condition, which is saved in `data/initial_condition.mat` and will be loaded automatically.
  The simulation results for the coordinates of the telomere loci and forces over time are saved in `data/yadd_*`
  User can use `VECTOR_package.m` to generate visualizations of the simulation results.

* The partial differential equations (PDEs) are solved using finite differencing. Time integration is done using the implicit method via Matlab's `ode15s` where linear equations are solved using the direct method. The jacobian is an estimate of the full jacobian based on the block diagonals of the  concentration field and the telomere degrees of freedom respectively while setting the off block diagonals to zero. Users are encouraged to use iterative (Krylov) method which requires modification of `ode15s` and can substantially speed up the simulations.