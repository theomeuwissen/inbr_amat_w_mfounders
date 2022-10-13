function ML1992(ped)
n=size(ped,1)+1   #ped(i,1:2)=[sire dam]; 0 for missing
point=zeros(Int,n)
F=zeros(n)
L=zeros(n)
D=zeros(n)
PED=zeros(Int,n,2)           #create space for 1-group
PED[2:n,1:2]=ped[1:n-1,1:2].+1  

F[1]=-1.   #actually F[0]
for i=2:n
  is=PED[i,1]
  id=PED[i,2]
  PED[i,1]=max(is,id)
  PED[i,2]=min(is,id)
  D[i]=.5-.25*(F[is]+F[id])
  if(is==1)&&(id==1)   #base animal
    F[i]=0.0
  else  #non-base
    Fi=-1.0
    L[i]=1.0
    j=i   #youngest ancestor of i in llist
    while (j>1) #loop through llist
      r=L[j]/2
      (ks,kd)=PED[j,1:2]
      k=j
      if(ks>1)  #include parents in ancestor list
        while (point[k]>ks)  #find slot in llist for ks
	  k=point[k]
	end
	L[ks]+=r
	if (ks != point[k])   #incl ks in llist
	  point[ks]=point[k]  #point[ks]=>next anc in llist
	  point[k]=ks         #point(previous anc)=>ks
	end
	if(kd>1)
	  while(point[k]>kd)  #find slot in llist for kd
	    k=point[k]
	  end
	  L[kd]+=r
	  if(kd !=point[k])     #incl kd in llist
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
    F[i]=Fi
  end #non-base animal
end #animal loop

return F[2:n],D[2:n]
end #function
