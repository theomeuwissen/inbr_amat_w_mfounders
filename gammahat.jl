# estimates relationshio matrix of metafounders GAMMA
# given a relationship matrix G of animals and a contribution matrix Q of metafounders to the animals in G
# estimated by the statistical model: y = X*b +e
# where y = data vector of the elements of G (lower triangle is used)
#       b = data vector of the estimates of GAMMA (lower triangle)
#       X_1 = row 1 of the coefficient matrix
#belonging to  GAM11   GAM21      GAM22   GAM31     GAM32     GAM33   GAM41    .... 
#       X_1   =[Q11^2  2*Q12*Q11  Q12^2  2*Q13*Q11  2*Q13*Q12 Q13^2   2*Q14*Q11 ...
#       X_2   =[Q11Q21 (Q12Q21+   Q12Q22 (Q13Q21+   (Q13Q22+  Q13Q23  (Q14Q21+
#                       Q22Q11)           Q23Q11)    Q23Q12)           Q24Q11)
#

function gammahat(G,Q)
N=size(Q,1)
Nmf=size(Q,2)
Neq=Int(Nmf*(Nmf+1)/2)
bhat=zeros(Neq)
GAMMA=zeros(Nmf,Nmf)
XX=zeros(Neq,Neq)
Xy=zeros(Neq)
X=zeros(10000,Neq)  #evaluate 1000 data-points at a time
y=zeros(10000)

# set up normal equations XX and Xy
l=0
for i=1:size(Q,1), j=1:i  #loop over data-points
    l+=1
    y[l]=G[i,j]
    ll=0
    for ii=1:Nmf, jj=1:ii  #loop over elements of GAMMA
       ll+=1
       if(ii==jj)
            X[l,ll]=Q[i,ii]*Q[j,jj]
       else
            X[l,ll]=Q[i,ii]*Q[j,jj]+Q[j,ii]*Q[i,jj]
       end    
    end
    if(l==10000)||(i==j==N)  #X is full or done
#      println(X[1:l,1:Neq])
#      println(y[1:l])
      XX+=X[1:l,1:Neq]'*X[1:l,1:Neq]
      Xy+=X[1:l,1:Neq]'y[1:l]
      l=0; 
    end
    if(i%10000==0)&&(i==j)
      println(" relationships of animal ",i)
    endif  
end

#solve equations
bhat=XX\Xy
l=0
for i=1:Nmf, j=1:i
  l+=1
  GAMMA[i,j]=bhat[l]
  GAMMA[j,i]=bhat[l]
end

return GAMMA
end #function
