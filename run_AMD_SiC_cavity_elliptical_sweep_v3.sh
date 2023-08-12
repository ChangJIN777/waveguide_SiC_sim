module load python/3.10.9-fasrc01 #Load python module
mamba activate wv_sol_env # activate the wavesolver environment
module load intel/23.0.0-fasrc01 
module load lumerical-seas/2021_7bf43e7149-fasrc01
python ./AMD_SiC_cavity_elliptical_optimization_v1.py
