# inbr_amat_w_mfounders

This repository contains some Julia functions to calculate inbreeding coeffients and/or the A-matrix.
Also Henderson's inverse of the A-matrix is calculated.
Calculations are performed with or without Metafounders (Legarra et al. 2015)

PED arrays consist of a sire and a dam columns of integers (2 columns).
Pedigrees are numbered 1,2,3,...,n and sorted from old to young.
If Nmf metafounders exist: they are numbered 1,..,Nmf and the real animals are numbered Nmf+1,...,n.

lst22 contains a list of animal IDs for which the A matrix is needed. (this is often called the A22 matrix, referring to the A matrix of genotyped animals in ssGBLUP). The lst22 list is assumed sorted from old to young animals. In A22 matrix output IDs have same order as in lst22.

GAMMA is matrix of (NmfxNmf) relationships between the metafounders. Assumed known beforehand.

## packages required:
Before running the routines use the LinearAlgebra package:
using LinearAlgebra, SparseArrays


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



```
a22=zeros(N22*(N22+1)/2)

fil=open("A22.mat")

a22=read!(fil,a22)

close(fil)
```


## A22mf.jl: calculates A-relationship matrix incl metafounders for list of animals (lst22)
Usage: F = A22(ped,lst22,GAMMA)
F is returned which is a vector of inbreeding coefficients (Meuwissen&Luo, 1992 algorithm).
Writes: A22.mat file which contains the lower triangle elements (double-precision) of the A22 matrix.
The elements can be read in Julia into a vector a22 by:


```
a22=zeros(N22*(N22+1)/2)

fil=open("A22.mat")

a22=read!(fil,a22)

close(fil)
```

## fastA22mf.jl: calculates A-relationship matrix incl metafounders for list of animals (lst22) faster than A22mf.jl
Usage: fastA22(ped,lst22,GAMMA,f)
f is (n-Nmf)x1 vector of inbreeding coefficients (input; may be calculated by ML1992MF.jl).
Writes: fastA22mf.mat file which contains the lower triangle elements (double-precision) of the A22 matrix.
The elements can be read in Julia into a vector a22 by:


```
a22=zeros(N22*(N22+1)/2)  #no_elements equals size_of_file fastA22.mat in bytes divided by 8  

fil=open("fastA22mf.mat")

a22=read!(fil,a22)

close(fil)
```



## fastA22.jl: calculates A-relationship matrix for list of animals (lst22) faster than A22.jl
Usage: fastA22(ped,lst22,f)
f is (nx1) vector of inbreeding coefficients (input; may be calculated by ML1992.jl).
Writes: fastA22.mat file which contains the lower triangle elements (double-precision) of the A22 matrix.
The elements can be read in Julia into a vector a22 by:


```
a22=zeros(N22*(N22+1)/2)  #no_elements equals size_of_file fastA22.mat in bytes divided by 8  

fil=open("fastA22.mat")

a22=read!(fil,a22)

close(fil)
```





## AI.jl: calculates Henderson's sparse inverse of A matrix
Usage: AI(ped,f)
f= (nx1) vector of inbreeding coefficients (input; use ML1992.jl to get f)
ped pedigree in 2 columns (sire dam; 0 for unknown)

Output: AI.mat is sparse lower triangle of A-inverse (3 columns: I J and AI(I,J))



## AI_mf.jl: calculates Henderson's sparse inverse of A matrix with metafouders
Usage: AI_mf(ped,f,GAMMA)
f= ((n-Nmf)x1) vector of inbreeding coefficients for the real animals only (input; use ML1992mf.jl to get f)

ped pedigree in 2 columns (sire dam). First Nmf entries are metafounders with '0 0' as parents. Real animals have IDs: Nmf+1,..,n

GAMMA= (NmfxNmf) relationship matrix of metafounders (input). Their inbreeding coefficients are diag(GAMMA)/2. 

Output: AI_mf.mat is sparse lower triangle of A-inverse including metafounder part (3 columns: I J and AI(I,J))




## gammahat.jl: estimates gamma matrix of metafounder relationships

Usage: gamma=gammhat(G,Q)

where: G is nxn (genomic) relationship-matrix or A22 matrix
Q = (nxNmf) matrix of genetic contribution of Nmf metafounders to n animals

Model used: G = Q * GAMMA * Q' + E  

The data-matrix G is vectorised and the unknowns (GAMMA) are also vectorised to yield a usual least squares model. 

Thus given a relationship matrix G (or A22) of animals and a contribution matrix Q of metafounders to the animals in G
the estimation is by the statistical model: y = X*b +e  

where y = data vector of the elements of G (lower triangle is used)  
       b = data vector of the estimates of GAMMA (lower triangle)  
       X_1 = row 1 of the coefficient matrix  
belonging to:
```
              GAM11   GAM21      GAM22   GAM31     GAM32     GAM33   GAM41    ....  
       X_1   =[Q11^2  2*Q12*Q11  Q12^2  2*Q13*Q11  2*Q13*Q12 Q13^2   2*Q14*Q11 ...  
       X_2   =[Q11Q21 (Q12Q21+   Q12Q22 (Q13Q21+   (Q13Q22+  Q13Q23  (Q14Q21+  
                       Q22Q11)           Q23Q11)    Q23Q12)           Q24Q11)  
```


