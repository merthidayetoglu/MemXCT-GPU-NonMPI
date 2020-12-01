#!/bin/bash

#DOMAIN INFORMATION
export NUMTHE=2400
export NUMRHO=2048
export PIXSIZE=1
#SOLVER DATA
export NUMITER=24
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
export THEFILE=./datasets/ADS4_theta.bin
export SINFILE=./datasets/ADS4_sinogram.bin
export OUTFILE=./datasets/recon_ADS4.bin

export OMP_NUM_THREADS=1

build/exe/src/main.cu.exe
