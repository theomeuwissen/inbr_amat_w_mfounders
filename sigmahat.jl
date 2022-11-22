# following Masuda et al. (2021) and references therein:
#  A* = A + Q*SIGMA*Q'
#  A* = ped relationships incl groups
# Q = group contribution matrix
# A = ped relationships without groups
# sigmahat estimates SIGMA given a general relationship matrix G (or A)
#
# i.e. the same model can be applied to the genomic relationship matrix G (we assume genomic relationships include groups resulting in SIGMA_G)
# 
# We can then extend Vitezica et al. (2011) to a model with groups:
# Gadjust = G + Q*(SIGMA_A - SIGMA_G)*Q' 
#
#
# thinning= 10 an integer indicating that every 10th element of G is to be used for the estimation of GAMMA (=0 means use all elements) 
#
# given a relationship matrix G of animals and a contribution matrix Q of metafounders to the animals in G
# estimated by the statistical model: y = X*b +e
# where y = data vector of the elements of G (lower triangle is used)
#       b = data vector of the estimates of SIGMA (lower triangle)
#       X_1 = row 1 of the coefficient matrix
#belonging to  SIGM11   SIGM21     SIGM22  SIGM31     SIGM32   SIGM33   SIGM41    .... 
#       X_1   =[Q11^2  2*Q12*Q11  Q12^2  2*Q13*Q11  2*Q13*Q12 Q13^2   2*Q14*Q11 ...
#       X_2   =[Q11Q21 (Q12Q21+   Q12Q22 (Q13Q21+   (Q13Q22+  Q13Q23  (Q14Q21+
#                       Q22Q11)           Q23Q11)    Q23Q12)           Q24Q11)
#
#
# 

function sigmahat(G,Q,thinning)
N=size(Q,1)
Nmf=size(Q,2)
Neq=Int(Nmf*(Nmf+1)/2)
bhat=zeros(Neq)
SIGMA=zeros(Nmf,Nmf)
XX=zeros(Neq,Neq)
Xy=zeros(Neq)
X=zeros(10000,Neq)  #evaluate 1000 data-points at a time
y=zeros(10000)
thinning=max(1,Int(thinning))

# set up normal equations XX and Xy
l=0; k=0
for i=1:size(Q,1)
iprint=1
for j=1:i  #loop over data-points
  k+=1
  if(k%thinning==0)
    l+=1
    y[l]=G[i,j]
    ll=0
    for ii=1:Nmf, jj=1:ii  #loop over elements of SIGMA
       ll+=1
       if(ii==jj)
            X[l,ll]=Q[i,ii]*Q[j,jj]
       else
            X[l,ll]=Q[i,ii]*Q[j,jj]+Q[j,ii]*Q[i,jj]
       end    
    end
    if(l==10000)   #||(i==j==N)  #X is full or done
#      println(X[1:l,1:Neq])
#      println(y[1:l])
      XX+=X[1:l,1:Neq]'*X[1:l,1:Neq]
      Xy+=X[1:l,1:Neq]'y[1:l]
      l=0; 
    end
    if(i%10000==0)&&(iprint==1)
      println(" relationships of animal ",i); iprint=0
    end
  end  #thinning
end
end
if(l>0)  #add last records
      XX+=X[1:l,1:Neq]'*X[1:l,1:Neq]
      Xy+=X[1:l,1:Neq]'y[1:l]
end


#solve equations
#bhat=XX\Xy   #often many singularities
# use eigen-decomposition of XX=>  E*V*E'*bhat=Xy => bhat=E*V^-1*E'*Xy (use only first k-components of decomposition)
# use SVD for eigendecomposition of symmetric matrix
(E,V,Et)=svd(XX)
k=count(V.>1.0e-6)
bhat=E[:,1:k] * diagm(1 ./V[1:k]) * Et[:,1:k]'

l=0
for i=1:Nmf, j=1:i
  l+=1
  SIGMA[i,j]=SIGMA[j,i]=bhat[l]
end

return XX,Xy,SIGMA
end #function
