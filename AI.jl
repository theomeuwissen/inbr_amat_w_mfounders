#Henderson's A inverse
#F are inbreeding coefficients (input)

#Sparse lower triangle is writen to file AI.mat (3 columns: i j AI[i,j])

function AI(ped,F)
n=size(ped,1)
AI=spzeros(n,n)

for i=1:n
  (is,id)=ped[i,1:2]
  if(is==0)&&(id==0)
    AI[i,i]=1.0
  elseif (is>0)&&(id>0)
    d=.5-.25*(F[is]+F[id])
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
    d=.75-.25*F[is]
    AI[i,i]=1.0/d
    AI[is,is]+=.25/d
    AI[i,is]-=.5/d
    AI[is,i]-=.5/d
  elseif(id>0)
    d=.75-.25*F[id]
    AI[i,i]=1.0/d
    AI[id,id]+=.25/d
    AI[i,id]-=.5/d
    AI[id,i]-=.5/d
  end #if
end #for
#return AI
fil=open("AI.mat","w")
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
