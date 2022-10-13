#A22 calculates A matrix from the 'lst22' list of individuals
#IDs in pedigree are assumed consecutative numbered from old to young; (0 for missing sire/dam)
#lst22 is assumed sorted from low to high
# writes lower traingle in binary stream to A22.mat (Float64) [length(lst22)*(length(lst22)+1)/2 figures]
# returns inbreeding coefficients
# Theo Meuwissen Oct2022

function A22(ped,lst22)
n=size(ped,1)+1   #ped(i,1:2)=[sire dam]; 0 for missing
point=zeros(Int,n)
F=zeros(n)
L=zeros(n)
D=zeros(n)
PED=zeros(Int,n,2)           #create space for 1-group
PED[2:n,1:2]=ped[1:n-1,1:2].+1  #increase ID numbers by 1
lst22=lst22.+1               #increase ID numbers by 1
ANC_A22=Vector{Vector{Int32}}(undef,length(lst22))
L_A22=Vector{Vector{Float32}}(undef,length(lst22))

# count how many ancestors typically in lst22 animals
j=lst22[end]
ncnt=1
while (j>1)
  k=j
  (is, id)=PED[k,1:2]
  ks=max(is,id)
  kd=min(is,id)
  if(ks>0)
    while (point[k]>ks)
      k=point[k]
    end
    if(ks!=point[k])
      point[ks]=point[k]
      point[k]=ks
    end
    if(kd>0)
      while(point[k]>kd)
        k=point[k]
      end
      if(kd!=point[k])
         point[kd]=point[k]
	 point[k]=kd
      end
    end #if
  end #if
  ncnt+=1
  k=j
  j=point[j]
  point[k]=0
end #while
println(" number of ancestors of ",lst22[end]," = ",ncnt)
## finished counting ancestors


for i=1:length(lst22)
  ANC_A22[i]=zeros(Int32,0)
  L_A22[i]=zeros(Float32,0)
  sizehint!(ANC_A22[i],ncnt)
    sizehint!(L_A22[i],ncnt)
end  

# use Meuwissen_Luo_1992 algorithm for F and D
# store ancestor lists and 
i22=0   #volgnr in lst22
F[1]=-1.   #actually F[0]
for i=2:n
  if(i in lst22)
    i22+=1; in22=1
  else
    in22=0   #not in lst22
  end  
  is=PED[i,1]
  id=PED[i,2]
  PED[i,1]=max(is,id)
  PED[i,2]=min(is,id)
  D[i]=.5-.25*(F[is]+F[id])
  if(is==1)&&(id==1)   #base animal
    F[i]=0.0; 
    if(in22 == 1)        
        append!(ANC_A22[i22],i) #store i in ancestor list of i22
        append!(L_A22[i22],1.0)  #store contribution also 
    end
  else  #non-base
    Fi=-1.0
    L[i]=1.0
    j=i   #youngest ancestor in llist
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
      k=j                  #remember old ancest
      j=point[j]           #next youngest ancest

      if(in22 == 1)
        append!(ANC_A22[i22],k) #store k in ancestor list of i22
        append!(L_A22[i22],L[k])  #store contribution also 
      end

      L[k]=0.0             #clear
      point[k]=0           #clear
    end #while
    F[i]=Fi; 
  end #non-base animal
end #animal loop

# calculate A22 relationships
fid=open("A22.mat","w")   #write lower triangle in binary stream (Float64)
for i=1:length(lst22)
  for j=1:i-1
    Aij=0.0
    nexti=iterate(ANC_A22[i])  #loop over ancestors of i and j
    nextj=iterate(ANC_A22[j])
    while (nexti!==nothing) && (nextj!==nothing)
      if (nexti[1]==nextj[1])
         Aij+=L_A22[i][nexti[2]-1]*L_A22[j][nextj[2]-1]*D[nexti[1]]
	 nexti=iterate(ANC_A22[i],nexti[2]);
	 nextj=iterate(ANC_A22[j],nextj[2])
      elseif (nexti[1]>nextj[1])  #take next youngest ANC in i list
         nexti=iterate(ANC_A22[i],nexti[2]); 
      elseif (nextj[1]>nexti[1])  #take next youngest ANC in j list
	 nextj=iterate(ANC_A22[j],nextj[2])
      end #if
    end #while
    write(fid,Aij)
  end #j 
  ii=ANC_A22[i][1]  #ID of i
  write(fid,(1.0+F[ii]))   #write diagonal
end #i
close(fid)

return F[2:n]
end #function
