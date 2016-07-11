#In umean_profiles_v2.jl
include("perturb.jl")
function umean_profiles_v2(yplus,uplus,Np)
#----yplus=wall normal positions without perturbate
#----Np= number of assembles Np=0, master profile
#----Output paramters
#----ustep=row vector of u velocities of 5481 elements
#----ystep=row vector of y+ wall normal positions of 5481 elements
#----Initialize vectors
ypro=zeros(length(yplus));
ywozero=zeros(length(yplus));
ypperturbated=zeros(length(yplus));
gp=5481;                              #--Grid Points
ustep=zeros(gp);
U=zeros(gp);
#--------------------------------
ypperturbated=perturb(yplus[:],Np);
ypro=round(Int,ypperturbated);
while length(unique(ypro))<length(ypro)     #---check for repeated elements
     ypperturbated=perturb(yplus[:],Np);
     ypro=round(Int,ypperturbated);
end                    
sorty,sortu=sortyv(ypro[:],uplus);
ywozero=sorty[:];
unshift!(sorty,0);                             #--add zero to this array like first element warning
Nvf=length(sorty); 
#print(" sorty, ywozero= ",sorty, ywozero);
#-------Fill the step profile
for i=1:Nvf-1                                  
    U[sorty[i]+1:sorty[i+1]]=sortu[i];
end
#print("\n")
#print(" uplus= ", sortu);
#print("sorty= ", sorty)
return U, ywozero, sortu
end
#-----Function to sort vortical fissure y+
function sortyv(yraw,uraw)
#----Input
#----yraw   raw vector y without sort  (pertubated yp)
#----uraw   raw vector u without sort
#----Output
#----ysort  sorted vector y in ascending order
#----usort  sorted vector u according with original indexes in yraw
index=sortperm(yraw);
ysort=sort(yraw);
usort=uraw[index];
# print("y+ sorted",ysort,index,usort)
return ysort,usort
end
