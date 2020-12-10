#!/bin/bash

#DOMAIN INFORMATION
export NUMTHE=1500
export NUMRHO=1024
export PIXSIZE=1
#SOLVER DATA
export NUMITER=25
#TILE SIZE (MUST BE POWER OF TWO)
export SPATSIZE=128
export SPECSIZE=128
#BLOCK SIZE
export PROJBLOCK=256
export BACKBLOCK=256
#BUFFER SIZE
export PROJBUFF=48 #KB
export BACKBUFF=48 #KB
#I/O FILES
export THEFILE=./datasets/ADS3_theta.bin
export SINFILE=./datasets/ADS3_sinogram.bin
export OUTFILE=./datasets/recon_ADS3.bin

export OMP_NUM_THREADS=1

export PINDEX=/mnt/nvme0/memxctdata/ADS3/pidxfile.bin
export PVALUE=/mnt/nvme0/memxctdata/ADS3/pvalfile.bin
export BINDEX=/mnt/nvme0/memxctdata/ADS3/bidxfile.bin
export BVALUE=/mnt/nvme0/memxctdata/ADS3/bvalfile.bin


#nvprof --analysis-metrics -f -o analysis.nvvp build/exe/src/main.cu.exe
build/exe/src/main.cu.exe 5
