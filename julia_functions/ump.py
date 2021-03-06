from per import perturb
import numpy as np
def umean_profiles(yplus,uplus,Np): #program does nothing as written
#----yplus=wall normal positions without perturbate
#----Np= number of assembles Np=0, master profile
#----Output paramters
#----ustep=row vector of u velocities of 5481 elements
#----ystep=row vector of y+ wall normal positions of 5481 elements
    #print("yp2=", yplus)
    ypro=np.zeros(len(yplus));
    ypperturbated=np.zeros(len(yplus));
    ypperturbated=perturb(yplus[:],Np);
    #print("yp3per=", ypperturbated)
    ypro=ypperturbated.astype(np.int64)[:];           #-convert from float to int64
    #-----Check for repeated elements
    while len(set(ypro))<len(ypro):
          ypperturbated=perturb(yplus[:],Np);
          ypro=ypperturbated.astype(np.int64)[:];
    #---------------------------------
    sorty,sortu=sortyv(ypro[:],uplus);
    #print("sorty,sortu",sorty,sortu)
    sorty=np.append(0,sorty);
    #print("sorty=", sorty)
    Nvf=len(sorty)-1;                     #-Index of last element in sorty
    gp=5481;                              #--Grid Points
    ustep=np.zeros(gp);
    U=np.zeros(gp);
    for i in range(0,Nvf):
       U[sorty[i]:sorty[i+1]]=sortu[i];  #it overwrites steps, fix this
    return U

def sortyv(yraw,uraw):
   #----Input
   #----yraw   raw vector y without sort  (pertubated yp)
   #----uraw   raw vector u without sort
   #----Output
   #----ysort  sorted vector y in ascending order
   #----usort  sorted vector u according with original indexes in yraw
   index=np.argsort(yraw);
   ysort=np.sort(yraw);
   usort=uraw[index];
  # print("y+ sorted",ysort,index,usort)
   return [ysort, usort]
