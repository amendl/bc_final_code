#! /bin/sh

#SBATCH -p q_fs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 0-10:00
#SBATCH -o std_methanol.out
#SBATCH -e std_methanol.err


/home/users/mendla/custom_software/env/bin/python from_dimers_CC.py
