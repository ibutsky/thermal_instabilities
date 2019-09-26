
Purpose: to generate initial conditions necessary to run a thermal instability box with Enzo

Directions: 

1. Edit the top section of constants_and_parameters.py to set gas parameters
2. Run $python generate_initial_conditions.py. This will generate the final
   parameter file that's read in by Enzo (ThermalInstability.enzo), as well 
   as two necessary input files (cool_rates.in and perturbation.in). 
3. That's it! The file "submit.sbatch" has an example of how to run the 
   simulation from the command line
