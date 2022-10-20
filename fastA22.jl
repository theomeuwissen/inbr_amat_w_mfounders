# fastA22 calculates A matrix from the 'lst22' list of individuals
#IDs in pedigree are assumed consecutative numbered from old to young; 0 => unknown
#lst22 is assumed sorted from low to high
# f = input Nx1 vector of inbreeding coefficients 
# writes lower triangle in binary stream to fastA22.mat (Float64) [length(lst22)*(length(lst22)+1)/2 figures]
# Theo Meuwissen Oct2022

function fastA22(PED,lst22,F)
n=size(PED,1)   #ped(i,1:2)=[sire dam]; 0 for missing
n22=size(lst22,1)
in_lst22=zeros(Int32,n)   #indicator whether animal is in list
in_lst22[lst22].=1
point=zeros(Int,n)
L=zeros(n)
D=zeros(n)



# use Meuwissen_Luo_1992 algorithm for F and D
# store ancestor lists and 

# set up D
for i=1:n
  is=PED[i,1]
  id=PED[i,2]
  PED[i,1]=max(is,id)
  PED[i,2]=min(is,id)
  if(is>0)&&(id>0)
    D[i]=.5-.25*F[is]-.25*F[id]   # .25*(facts*(1-F[is])+factd*(1-F[id]))
  elseif(is>0)
    D[i]=.75-.25*F[is]
  elseif(id>0)
    D[i]=.75-.25*F[id]
  else
    D[i]=1.0
  end  
end
#println("D= ",D)

## calculate A22 relationships
#calculate column i of A22 as A*Ei where Ei is unit vector (of zeros with a one at position i)
#A*Ei=(I-P)^{-1}*D*(I-P')^{-1} * Ei
#where P-matrix has rows of i with elements-j equal 0.5 when j is parent and L=(I-P)^{-1}
#STEP 1: calculate column i of (I-P')^{-1} = row of L
#STEP 2: calculade D*(rowi_of_L)
#STEP 3 solve (I-P)*sol=L (store sol in L)



fid=open("fastA22.mat","w")   #write lower triangle in binary stream (Float64)
for ii=1:n22
    i=lst22[ii]
# step 1,2: get row of L and multiply with D
    L=zeros(n); 
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
      L[k]*=D[k]
      point[k]=0           #clear
    end #while

#step 3 solve (I-P)*sol=L (store sol in L)
    for k=1:i  #founders dont have pedigree; dont need relationships beyond i
       (ks,kd)=PED[k,1:2]
       if(ks>0)
        L[k]+=.5*L[ks]
       end
       if(kd>0)
        L[k]+=.5*L[kd]
       end
       if(in_lst22[k]>0)
          write(fid,L[k])
       end    	
     end #k
     println("individ ",lst22[ii])
end #lst22 loop
close(fid)


return 
end #function
