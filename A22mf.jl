# A22mf calculates A matrix from the 'lst22' list of individuals
#IDs in pedigree are assumed consecutative numbered from old to young; (first size(GAMMA) id-numbers are metafounders)
#lst22 is assumed sorted from low to high
# writes lower traingle in binary stream to A22mf.mat (Float64) [length(lst22)*(length(lst22)+1)/2 figures]
# returns inbreeding coefficients (except for metafounders)
# Theo Meuwissen Oct2022

function A22mf(PED,lst22,GAMMA)
n=size(PED,1)   #ped(i,1:2)=[sire dam]; 0 for missing
Nmf=size(GAMMA,1) # no. of metafounders
point=zeros(Int,n)
F=zeros(n)
L=zeros(n)
D=zeros(n)
ANC_A22=Vector{Vector{Int32}}(undef,length(lst22))
L_A22=Vector{Vector{Float32}}(undef,length(lst22))

# count how many ancestors typically in lst22 animals
j=lst22[end]
ncnt=1
while (j>0)
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
println(" Nmetafounders= ",Nmf)

for i=1:length(lst22)
  ANC_A22[i]=zeros(Int32,0)
  L_A22[i]=zeros(Float32,0)
  sizehint!(ANC_A22[i],ncnt)
  sizehint!(L_A22[i],ncnt)
end  

# use Meuwissen_Luo_1992 algorithm for F and D
# store ancestor lists and 
for i=1:Nmf
  F[i]=.5*GAMMA[i,i]
end  
i22=0   #volgnr in lst22
for i=Nmf+1:n
  if(i in lst22)
    i22+=1; in22=1
  else
    in22=0   #not in lst22
  end  
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
    j=i   #youngest ancestor in llist
    while (j>0) #loop through llist
      r=L[j]/2
      (ks,kd)=PED[j,1:2]
      k=j
      if(ks>0)  #include parents in ancestor list
        while (point[k]>ks)  #find slot in llist for ks
	  k=point[k]
	end
	L[ks]+=r
	if (ks != point[k])   #incl ks in llist
	  point[ks]=point[k]  #point[ks]=>next anc in llist
	  point[k]=ks         #point(previous anc)=>ks
	end
	if(kd>0)
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
      k=j                  #remember old ancest
      j=point[j]           #next youngest ancest

      if(in22 == 1)
        append!(ANC_A22[i22],k) #store k in ancestor list of i22
        append!(L_A22[i22],L[k])  #store contribution also 
      end
      if(k>Nmf)
         Fi+=D[k]*L[k]*L[k]
	 L[k]=0.0     #clear
      else
         Fi+=L[k]*GAMMA[k,1:Nmf]'*L[1:Nmf]
      end	 
      point[k]=0           #clear
    end #while
    L[1:Nmf]=zeros(Nmf)    
    F[i]=Fi; 
end #animal loop

# calculate A22 relationships
fid=open("A22mf.mat","w")   #write lower triangle in binary stream (Float64)
for i=1:length(lst22)
  for j=1:i-1
    Aij=0.0; Li=zeros(Nmf); Lj=zeros(Nmf)
    nexti=iterate(ANC_A22[i])  #loop over ancestors of i and j
    nextj=iterate(ANC_A22[j])
    while (nexti!==nothing) || (nextj!==nothing)
      if(nexti!==nothing) && (nextj!==nothing) #outer
       if(nexti[1]>Nmf) && (nextj[1]>Nmf)  #inner
        if (nexti[1]==nextj[1])
         Aij+=L_A22[i][nexti[2]-1]*L_A22[j][nextj[2]-1]*D[nexti[1]]
	 nexti=iterate(ANC_A22[i],nexti[2]);
	 nextj=iterate(ANC_A22[j],nextj[2])
        elseif (nexti[1]>nextj[1])  #take next youngest ANC in i list
         nexti=iterate(ANC_A22[i],nexti[2]); 
        elseif (nextj[1]>nexti[1])  #take next youngest ANC in j list
	 nextj=iterate(ANC_A22[j],nextj[2])
        end #if
       else #inner
        if(nexti[1]<=Nmf)
          Li[nexti[1]]=L_A22[i][nexti[2]-1]
	end  
        if(nextj[1]<=Nmf) 
	  Lj[nextj[1]]=L_A22[j][nextj[2]-1]
	end  
 	nexti=iterate(ANC_A22[i],nexti[2]);
	nextj=iterate(ANC_A22[j],nextj[2])
       end #inner	
     elseif(nexti!==nothing) #outer
        if(nexti[1]<=Nmf)
          Li[nexti[1]]=L_A22[i][nexti[2]-1]
	end  
 	nexti=iterate(ANC_A22[i],nexti[2]);
     elseif(nextj!==nothing) #outer       
        if(nextj[1]<=Nmf) 
	  Lj[nextj[1]]=L_A22[j][nextj[2]-1]
	end  
	nextj=iterate(ANC_A22[j],nextj[2])
     end #outer
    end #while
    Aij+=sum(Li'*GAMMA*Lj)  #contribution of metafounders
    write(fid,Aij)
  end #j
  nexti=iterate(ANC_A22[i])  #find ID of i
  ii=nexti[1]  #ID of i
  println(" F(i)= ",ii," ",F[ii])
  write(fid,(1.0+F[ii]))   #write diagonal
end #i
close(fid)

return F[Nmf+1:n]
end #function
