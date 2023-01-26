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


## renum_ped.jl: renumbers pedigree
 usage: (ped, irenum) = renum_ped(PED,select)   
input:   
 PED (input) contains 3 columns: ID SIRE DAM (sorted from old to young, i.e. sires/dams come before offspring)   
 SIRE == 0 or not in ID column: implies missing SIRE (if SIRE not in ID column but not 0: SIREid is treated as group)   
 select: (vector) selects set of animals from pedigree (includes their ancestors); select=[] selects all IDs  
   
output:   
 ped (output) = 2 column integer matrix of renumbered IDs (sire and dam) (new IDs equal row number and are thus omitted)    
                in case of groups: first Ngroup IDs belong to groups and have "0 0" parents; followed by animalIDs who may have parent/group IDs (>0)   
 irenum = vector of original IDs (row-numbers refer to original IDs)   
                in case of groups: the first Ngroup IDs are groups    





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




## sigmahat.jl: estimates sigma matrix of Masuda (2021)

Usage: sigma=sigmahat(G,Q,thinning)   

 following Masuda et al. (2021) and references therein:  
  A* = A + Q*SIGMA*Q'   
  A* = ped relationships incl groups (usually A22 matrix)   
 Q = (nxNmf) group contribution matrix   
 A = ped relationships without groups (usually A22 matrix)      
 sigmahat estimates SIGMA given a general relationship matrix G (or A)    

I.e. the same model can be applied to the genomic relationship matrix G (we assume genomic relationships include groups resulting in SIGMA_G)
    
 We can then extend Vitezica et al. (2011) to a model with groups:   
 Gadjust = G + Q*(SIGMA_A - SIGMA_G)*Q'   
    
    
thinning= 10: an integer indicating that every 10th element of G is to be used for the estimation of GAMMA (saves time)    
        =0 means use all elements  of G    
    
    
   

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



## upg_relations.jl : estimates relationships between unknown parent groups (GAMMA)
 usage: (GAMMA,F)=upg_relations(PED3,Nupg)   
 Nupg=number of UPGs  
 PED3(1:n,1:3)  contains:  
 for row 1,..,Nupg : all (integer) elements are 0  
 for the real animals row Nupg+1,..,n: a renumbered pedigree with row i:  SIREi DAMi UPGi  
 if SIREi is unknown use UPG_of_sire; same for DAMi; if UPGi is unknown: use 0   
 pedigree is renumbered consecutivel for animals: Nupg+1,Nupg+2,...,n    (renumbered from old to young animals)   
 NOTE: UPGs are also renumbered consecutively: 1,..,Nupg   
 Estimation of diagonals of GAMMA: 2*F(animals_in_same_UPG)   
 Offdiagonals of GAMMA: GAMMA(i,j)=min(GAMMA(i,i),GAMMA(j,j))   (following Aguilar & Misztal, 2008).   
 uses: include("ML1992MF.jl") (to calculate F); using Statistics

##  Atimesx.jl: calculates A*x without setting up A (algorithm by Colleau et al (2002))
 usage s=Atimesx(x,ped,F)   
 x=vector of reals   
 ped = pedigree (2 column: column of sires; column of dams); 0 = missing parents   
 F is vector of inbreeding coefficients   
   unknown Ax = A*x= (I-L)^-1 *D* (I-L')^-1 *x    

 output : vector of A*x


## impute_parents: randomly imputes missing parents from their UPG group
usage PED = impute_parents(ped)    

 ped = 3 column pedigree: sire dam UPG    
 UPG = Unknown Parent Group = coded as negative integer ID-numbers (or 0 if unknown)   
 Unknown parents in columns 1 and 2 can be identified by their UPG-id (or 0 if true founder animals)   
 impute_parents(): for any missing parent randomly select a parent from the same UPG   
 Note: by having different UPGs for males/sires and females/dams we can avoid that a male is used as dam or vice versa.   

 output: PED is pedigree where true unknown parents (indicated 0) remain unknown but all other parents are randomly imputed from their UPG


