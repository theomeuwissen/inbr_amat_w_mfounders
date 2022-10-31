# upg_relations.jl : estimates relationships between unknown parent groups (GAMMA)
# usage: (GAMMA,F)=upg_relations(PED3,Nupg)
# Nupg=number of UPGs 
# PED3(1:n,1:3)  contains:
# for row 1,..,Nupg : all (integer) elements are 0 
# for the real animals row Nupg+1,..,n: a renumbered pedigree with row i:  SIREi DAMi UPGi  
# if SIREi is unknown use UPG_of_sire; same for DAMi; if UPGi is unknown: use 0
# pedigree is renumbered consecutivel for animals: Nupg+1,Nupg+2,...,n    (renumbered from old to young animals)
# NOTE: UPGs are also renumbered consecutively: 1,..,Nupg 
# Estimation of diagonals of GAMMA: 2*F(animals_in_same_UPG)
# Offdiagonals of GAMMA: GAMMA(i,j)=min(GAMMA(i,i),GAMMA(j,j))   (following Aguilar & Misztal, 2008).
# uses: include("ML1992MF.jl") (to calculate F); using Statistics


function upg_relations(PED3,Nupg)
#using Statistics
include("ML1992MF.jl")
n=size(PED3,1)
GAMMA=zeros(Nupg,Nupg)

F=zeros(n-Nupg)
minbold=0.0
for iter=1:100
  F=ML1992MF(PED3[1:n,1:2],GAMMA)
  avginb=zeros(Nupg)
  navginb=zeros(Int64,Nupg)
  for i=Nupg+1:n
    if(PED3[i,1]>Nupg)&&(PED3[i,2]>Nupg)&&(PED3[i,3]>0)
       avginb[PED3[i,3]]+=2*F[i-Nupg]
       navginb[PED3[i,3]]+=1
    end
  end
  minb=sum(F)/size(F,1)
  for i=1:Nupg
    if(navginb[i]>0)
      GAMMA[i,i]=avginb[i]/navginb[i]
    end
    for j=1:i-1
      GAMMA[i,j]=min(GAMMA[i,i],GAMMA[j,j])
      GAMMA[j,i]=GAMMA[i,j]
    end
  end
  println(" iteration_criterion= ",abs(minbold-minb))
  if(abs(minbold-minb)<1.E-6)
     break
  end
  minbold=minb
end
return GAMMA,F
end #function