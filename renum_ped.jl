# renum_ped.jl renumbers pedigree
# usage: (ped,irenum) = renum_ped(PED,list)
#input:
# PED (input) contains 3 columns: ID SIRE DAM (sorted from old to young, i.e. sires/dams come before offspring)
# SIRE == 0 or not in ID column: implies missing SIRE (if SIRE not in ID column but not 0: SIREid is treated as group)
# select =  vector of IDs that are to be included in the pedigree (their ancestors will be included); use "select=[]" to include all IDs 
#output:
# ped (output) = 2 column integer matrix of renumbered IDs (sire and dam) (new IDs equal row number and are thus omitted)
#                in case of groups: first Ngroup IDs belong to groups and have "0 0" parents; followed by animalIDs who may have parent/group IDs (>0)
# irenum = vector of original IDs (row-numbers refer to original IDs)
#                in case of groups: the first Ngroup IDs are groups   


function renum_ped(PED,select)
N=size(PED,1)

# add ancestors to select list
if(length(select)>0)
  for i=N:-1:1
     if (PED[i,1] in select)
       if !(PED[i,2] in select); push!(select,PED[i,2]); end #add to select
       if !(PED[i,3] in select); push!(select,PED[i,3]); end #add to select
     end
  end     
end

# renumber pedigree
list=Dict{Any, Int}() #empty dictionary
groups=Dict{Int,Any}() #opposite keypairing to find group-IDs back
n=0  #number of animals
m=0  #number of groups (negative number and starts at 0)
for i=1:N
 if(length(select)==0)||(PED[i,1] in select)
  if !(haskey(list,PED[i,2]))
     list[PED[i,2]]=m
     groups[m]=PED[i,2];  m-=1
  end
  if !(haskey(list,PED[i,3]))
     list[PED[i,3]]=m;  
     groups[m]=PED[i,3];  m-=1
  end   
  n+=1
  list[PED[i,1]]=n
 end 
end


ntot=length(list)
if(m==-1) #actually no groups
  m=0
  ntot=ntot-1
end
ped=zeros(Int,ntot,2)
irenum=Vector{Any}(undef,ntot)
ii=0
for i=m+1:N
  ii+=1
  if(i>0)&&(length(select)==0 || PED[i,1] in select)
    ped[ii,1]=list[PED[i,2]]-m #note: m<=0
    ped[ii,2]=list[PED[i,3]]-m
    irenum[ii]=PED[i,1]
  elseif(i<=0)
    irenum[ii]=groups[i]
  else
    ii-=1
  end  
end  


return ped,irenum

end #function