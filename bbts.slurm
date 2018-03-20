#!/bin/bash 
#SBATCH -t 1:00:00 
#SBATCH -N 1 
#SBATCH -n 16
time ./bbts -hist all -out all -g all -b azimuth -b offset -d vegas
