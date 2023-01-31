# PED = impute_parents(ped) 
#
# ped = 3 column pedigree: sire dam UPG
# UPG = Unknown Parent Group = coded as negative integer ID-numbers (or 0 if unknown)
# Unknown parents in columns 1 and 2 can be identified by their UPG-id (or 0 if true founder animals)
# impute_parents(): for any missing parent randomly select a parent from the same UPG
# Note: by having different UPGs for males/sires and females/dams we can avoid that a male is used as dam or vice versa.
#
# output: PED is pedigree where true unknown parents (indicated 0) remain unknown but all other parents are randomly imputed from their UPG


function impute_parents(ped)

PED=deepcopy(ped)

#set up lists of members UPG
Nupg=-findmin(ped[:,3])[1]
upg_members=Vector{Vector{Int64}}(undef,Nupg)  #Vector of vectors
for i=1:Nupg
  upg_members[i]=Int[]
end

# populate the member-list
for i=1:size(ped,1)
   iupg=-ped[i,3]
   if(iupg>0)   #known upg
      push!(upg_members[iupg],i)
   end
end

#impute unknown parents => random parents
for i=1:size(ped,1)
 for j=1:2 
  if(ped[i,j]<0)
    iupg=-ped[i,j]
    n=length(upg_members[iupg])
    if(n>0)
      is=upg_members[iupg][rand(1:n)]
      if(is==ped[i,j%2+1])  #sire==dam: sample again
        is=upg_members[iupg][rand(1:n)]
      end	
    else
      is=0
    end
    PED[i,j]=is
  end
 end
end

return PED
end #function