#! /bin/sh

#SBATCH -p q_fs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 0-10:00
#SBATCH -o std_NH3.out
#SBATCH -e std_NH3.err


/home/users/mendla/custom_software/env/bin/python from_dimers_NH3_CC.py
