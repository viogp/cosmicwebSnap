1. Generate input files for CUTE in *generate_input*, using the adequate model:
  - **cute_input_r.py** 
  - **cute_input_gflux_z.py**

2. Copy the input file to the *run_cute* directory:
   '''> cp [file] ../run_cute/in.txt'''

3. Run CUTE from the *run_cute* directory
   '''> ./qsub_cute.csh  '''

   with the adequate model and z/r-space option.

   This calls run_cute.csh and generates the adequate inpute parameter file within the *params* directory.

4. Plot
  - **xi_bias.py** : xi(r) and bias for each survey at a single z
  - **z_xi_bias.py** : xi(r) and bisa for different redshifts
  - **gflux_z_data.py** :  xi(s) compared to the Favole data.