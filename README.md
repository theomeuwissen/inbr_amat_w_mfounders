# inbr_amat_w_mfounders

This repository contains some Julia functions to calculate inbreeding coeffients and/or the A-matrix.
Also Henderson's inverse of the A-matrix is calculated.
Calculations are performed with or without Metafounders (Legarra et al. 2015)

PED arrays consist of a sire and a dam columns of integers (2 columns).
Pedigrees are numbered 1,2,3,...,n and sorted from old to young.
If Nmf metafounders exist: they are numbered 1,..,Nmf and the real animals are numbered Nmf+1,...,n.

lst22 contains a list of animal IDs for which the A matrix is needed. (this is often called the A22 matrix, referring to the A matrix of genotyped animals in ssGBLUP). The lst22 list is assumed sorted from old to young animals.

GAMMA is matrix of (NmfxNmf) relationships between the metafounders.

## packages required:
Before running the routines use the LinearAlgebra package:
using LinearAlgebra


## ML1992.jl: calculates Inbreeding coefficients.
Usage:  F = ML1992(ped)  
i.e. F is returned which is a vector of inbreeding coefficients
(Meuwissen&Luo, 1992 algorithm)


## ML1992MF.jl: calculates Inbreeding coefficients with metafounders.
Usage:   F = ML1992MF(ped,GAMMA)
i.e. F is returned which is a vector of inbreeding coefficients of real aniamls.
The relationship matrix of metafounders (GAMMA) is input and its size defines the numumber of metafounders.
(combined Meuwissen&Luo, 1992; and Vanraden, 1992; Legarra et al. 2015 algorithm)

## A22.jl: calculates A-relationship matrix for list of animals (lst22)
Usage: F = A22(ped,lst22)
F is returned which is a vector of inbreeding coefficients (Meuwissen&Luo, 1992 algorithm).
Writes: A22.mat file which contains the lower triangle elements (double-precision) of the A22 matrix.
The elements can be read in Julia into a vector a22 by:

'''
a22=zeros(N22*(N22+1)/2)
fil=open("A22.mat")
a22=read!(fil,a22)
close(fil)
'''


## A22mf.jl: calculates A-relationship matrix incl metafounders for list of animals (lst22)
Usage: F = A22(ped,lst22,GAMMA)
F is returned which is a vector of inbreeding coefficients (Meuwissen&Luo, 1992 algorithm).
Writes: A22.mat file which contains the lower triangle elements (double-precision) of the A22 matrix.
The elements can be read in Julia into a vector a22 by:

'''
a22=zeros(N22*(N22+1)/2)
fil=open("A22.mat")
a22=read!(fil,a22)
close(fil)
'''



