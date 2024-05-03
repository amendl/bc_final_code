#! /bin/sh

#SBATCH -p q_fs
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=30G
#SBATCH -t 0-10:00
#SBATCH -o std_Me.out
#SBATCH -e std_Me.err


/home/users/mendla/custom_software/env/bin/python from_dimers_Me_CC.py
