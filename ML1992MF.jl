# MeuwisenLuo1992 algorithm with MetaFounders (Legarra, 2015)
# ped: two columns: sires, dams; first Nmf entries are Metafounders (real animals start at position Nmf+1,..,N)
# GAMMA: (NmfxNmf) matrix of relationships between Metafounders

function ML1992MF(PED,GAMMA)
n=size(PED,1)   #ped(i,1:2)=[sire dam]; Metafounders have [0 0] as sires and dams
Nmf=size(GAMMA,1)
point=zeros(Int,n)
F=zeros(n)
L=zeros(n)
D=zeros(n)

for i=1:Nmf
  F[i]=.5*GAMMA[i,i]
end  


for i=Nmf+1:n
  is=PED[i,1]
  id=PED[i,2]
  PED[i,1]=max(is,id)
  PED[i,2]=min(is,id)
  facts=1; factd=1
  if(is<=Nmf); facts=2; end
  if(id<=Nmf); factd=2; end
  D[i]=.25*(facts*(1-F[is])+factd*(1-F[id])) #VanRaden1992
    Fi=-1.0
    L[i]=1.0
    j=i   #youngest ancestor of i in llist
    while (j>Nmf) #loop through llist
      r=L[j]/2
      (ks,kd)=PED[j,1:2]
      k=j
      if(ks>0)  #include parents in ancestor list
        while (point[k]>ks)  #find slot in llist for ks
	  k=point[k]
	end
	L[ks]+=r
	if (ks != point[k])&&(ks>Nmf)   #incl ks in llist
	  point[ks]=point[k]  #point[ks]=>next anc in llist
	  point[k]=ks         #point(previous anc)=>ks
	end
	if(kd>0)
	  while(point[k]>kd)  #find slot in llist for kd
	    k=point[k]
	  end
	  L[kd]+=r
	  if(kd !=point[k])&&(kd>Nmf)     #incl kd in llist
	    point[kd]=point[k]  #point[kd]=>next anc in llist
	    point[k]=kd         #point(previous anc)=>kd
	  end
	end
      end
      Fi+=D[j]*L[j]*L[j]
      L[j]=0.0             #clear
      k=j                  #store old ancest
      j=point[j]           #next youngest ancest
      point[k]=0           #clear
    end #while
    Fi+=sum(L[1:Nmf]'*GAMMA*L[1:Nmf]); L[1:Nmf]=zeros(Nmf)
    F[i]=Fi
end #animal loop

return F[Nmf+1:n]
end #function
