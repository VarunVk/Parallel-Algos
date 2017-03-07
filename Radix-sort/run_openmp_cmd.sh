#!/bin/sh
set -x

./rs_openmp /export/scratch/CSCI-5451/assignment-2/1M.txt 1 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/1M.txt 2 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/1M.txt 4 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/1M.txt 8 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/1M.txt 16

./rs_openmp /export/scratch/CSCI-5451/assignment-2/10M.txt 1 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/10M.txt 2 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/10M.txt 4 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/10M.txt 8 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/10M.txt 16 

./rs_openmp /export/scratch/CSCI-5451/assignment-2/100M.txt 1 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/100M.txt 2 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/100M.txt 4 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/100M.txt 8 
./rs_openmp /export/scratch/CSCI-5451/assignment-2/100M.txt 16 
