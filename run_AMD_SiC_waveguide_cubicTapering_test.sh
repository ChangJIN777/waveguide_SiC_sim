module load Anaconda3 #Load Anaconda module
source activate wv_sol_env # activate the wavesolver environment
module load intel/21.2.0-fasrc01
module load openmpi/4.1.1-fasrc01
module load lumerical-seas/2021_7bf43e7149-fasrc01
python ./AMD_SiC_waveguide_cubicTapering_test.py
