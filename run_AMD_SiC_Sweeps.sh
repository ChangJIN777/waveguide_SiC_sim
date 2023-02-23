#!/bin/bash
#SBATCH -c 10                # Number of cores (-c)
#SBATCH -t 0-00:30          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p serial_requeue   # Partition to submit to
#SBATCH --mem=100           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o myoutput_%j.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e myerrors_%j.err  # File to which STDERR will be written, %j inserts jobid
module load Anaconda3/5.3.0 #Load Anaconda module
module load lumerical-seas/2021_7bf43e7149-fasrc01
python3 AMD_SiNB_Sweeps.py
