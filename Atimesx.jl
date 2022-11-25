#  Atimesx.jl calculates A*x without setting up A
# usage s=Atimesx(x,ped,F)
# x=vector of reals
# ped = pedigree (2 column: column of sires; column of dams); 0 = missing parents
# F is vector of inbreeding coefficients
#   unknown s = A*x= (I-L)^-1 *D* (I-L')^-1 *x
#   step 1: solve for y in (I-L)'*y=x
#   step 2: z = D*y
#   step 3: solve for s in (I-L)*s=z
#   note1: if F(=inbreeding coefficients) present, D (=within fam var) is calculated (input D may be 0)
# output : vector of A*x

function Atimesx(x,ped,F)   

  n=size(x,1)

#  ! step 1: solve for y in (I-L)'*y=x
  y=deepcopy(x)
  d=zeros(n)      
  for i=n:-1:1
     if(ped[i,1]>0) &&  (ped[i,2]>0)
        d[i]=.5-.25*(F[ped[i,1]]+F[ped[i,2]])
     elseif(ped[i,1]>0)
        d[i]=.75-.25*F[ped[i,1]]
     elseif(ped[i,2]>0)
        d[i]=.75-.25*F[ped[i,2]]
     else
        d[i]=1.
     end 
     if(y[i]!=0.0)
        if(ped[i,1]>0); y[ped[i,1]]+=0.5*y[i]; end
        if(ped[i,2]>0); y[ped[i,2]]+=0.5*y[i]; end
     end 
  end

#  ! step 2: z = D*y
  y[1:n].*=d[1:n]
#  ! step 3: solve for s in (I-L)*s=z
  for i=1:n
     if(ped[i,1]>0)
        y[i]+=0.5*y[ped[i,1]]
     end 
     if(ped[i,2]>0)
        y[i]+=0.5*y[ped[i,2]]
     end 
  end
  
     return y
end #function Ax

