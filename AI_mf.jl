#Henderson's A inverse with metafounders (following VanRaden 1992; Legarra et al. 2015)
#F vector inbreeding coefficients (input; only for real animals)
#ped has two colmuns (sires dams). IDs 1..Nmf are Metafounders are included in ped and have (0 0) as parents
#    real animals start at Nmf+1,...,n
#
# GAMMA = (Nmf x Nmf) relationship matrix of metafounders

#Sparse lower triangle is writen to file AI.mat (3 columns: i j AI[i,j])

function AI_mf(ped,f,GAMMA)
n=size(ped,1)
Nmf=size(GAMMA,1)
F=[diag(GAMMA)/2.0    #full set of inbreeding coefficients
   f]
GAMinv=inv(GAMMA)

AI=spzeros(n,n)
AI[1:Nmf,1:Nmf]=GAMinv

for i=Nmf+1:n
  (is,id)=ped[i,1:2]
  if(is==0)&&(id==0)
    AI[i,i]=1.0
  elseif (is>0)&&(id>0)
    facts=factd=1.;
    if(is<=Nmf); facts=2.; end
    if(id<=Nmf); factd=2.; end
    d=.25*(facts*(1.0-F[is])+factd*(1.0-F[id]))
    AI[i,i]=1.0/d
    AI[is,is]+=.25/d
    AI[is,id]+=.25/d
    AI[id,is]+=.25/d
    AI[id,id]+=.25/d
    AI[i,is]-=.5/d
    AI[i,id]-=.5/d
    AI[is,i]-=.5/d
    AI[id,i]-=.5/d
  elseif(is>0)
    println("ERR: unknown dam should not occur")
  elseif(id>0)
    println("ERR: unknown sire should not occur")
  end #if
end #for
#return AI
fil=open("AI_mf.mat","w")
for i=1:n
  elems=findnz(AI[:,i])
  for j=1:size(elems[1],1)
     if(elems[1][j]<=i)
       write(fil,"$i $(elems[1][j]) $(elems[2][j]) \n")
     end
  end
end
#for i=1:n
#  for j=1:i
#     if(AI[i,j]!=0.0)
#     write(fil,"$i $j $(AI[i,j]) \n")
#     end
#  end
#end
close(fil)

end #function
