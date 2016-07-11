#In umean_profiles.jl
include("perturb.jl")
#include("perturbv2.jl")
include("perturb_log.jl")
include("interline.jl")
function umean_profiles(yplus::Vector,uplus::Vector,Np::Integer,gp::Integer,wvf::Int)
#----Input parameters  -------#
#----yplus=wall normal positions without perturbate
#----uplus=Mean velocities without pertrubate
#----Np= number of assembles Np=0, master profile
#----gp= number of grid points
#----wvf= width of vortical fissure
###---Output paramters  -------#                   
#----ustep=row vector of u velocities of 5481 elements
#----ystep=row vector of y+ wall normal positions of 5481 elements
#----Initialize vectors
const dy=wvf;
ypro=zeros(length(yplus));
ypperturbated=zeros(length(yplus));
#gp=5481;                              #--Grid Points
U=zeros(gp);    
#--------------------------------
ypperturbated=perturb(yplus[:],Np);
#ypperturbated=perturb_log(yplus[:],Np);        #----perturbation in log space
#ypperturbated=perturbv2(yplus[:],Np);         #--Negative and positive discriminator
ypro=round(Int,ypperturbated);
#print("ypro_1=",ypro)
#print("\n");
while length(unique(ypro))<length(ypro)     #---check for repeated elements
     ypperturbated=perturb(yplus[:],Np);
     #ypperturbated=perturb_log(yplus[:],Np);        #--pertrubation in log space
     #ypperturbated=perturbv2(yplus[:],Np);         #--Negative and positive discriminator
     ypro=round(Int,ypperturbated);
     #print("ypro_2=",ypro)
     #print("\n");
end                    
sorty,sortu=sortyv(ypro[:],uplus);
#print("sorty=",sorty)
#print("\n");
#print("uplus", uplus)
unshift!(sorty,0);                             #--add zero tho this array warning
#print("sorty=",sorty)
Nvf=length(sorty); 
#-------Fill the step profile
U[sorty[1]+1:sorty[2]]=sortu[1];              #--fill just first step
for i=2:Nvf-1                                 #--fill from 2nd step to last one
    if (sorty[i+1] > gp) 
      print("Vortical fissure y^+= ", sorty[i+1]);   
      print(" is out of the grid") 
      print("\n");
    end
    if dy>1                                  #--vortical fissures with different thickness
       U[sorty[i]+1:sorty[i]+dy]=interline(sortu[i-1],sortu[i],dy);
       print("y+= ", sorty[i]);
       print("\n");
    end     
    print("y+1= ", sorty[i+1]);
    print("\n");              
    U[sorty[i]+1+dy:sorty[i+1]]=sortu[i];
end
#print("sorty= ", sorty)
return U
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
