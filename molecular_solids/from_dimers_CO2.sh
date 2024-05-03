#! /bin/sh

#SBATCH -p q_fs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 0-10:00
#SBATCH -o std_CO2.out
#SBATCH -e std_CO2.err

/home/users/mendla/custom_software/env/bin/python from_dimers_CO2_2.py
