# initialize global variables used in your code
import numpy as np

def perturb(vectorr,Np):
    vector=vectorr[:];                                            #---copying by value
    if Np>0:
       Nve=len(vector);                                           #--vector length
       lindex=Nve-1;                                              #--Index last element of the vector, pyhton start at zero
       pert=0.4;
       h=np.zeros(Nve-1);
      # print("vector[lindex]=",vector[lindex])
      #print("Nve,vector",Nve,vector)
       for j in range(1,lindex):                                   #--perturb last and first position violate bc
       #  print("j=",j)
         h[j]=vector[j]-vector[j-1];                              #--Actual height of the y(j) step
         vector[j]=vector[j]+(pert*np.random.randn())*h[j];      #--gauss distr mu=0 sigma=1 (-5,5) normalized by 5 (check this for problems in b.c)
         while ((vector[j]<=vector[0]) | (vector[j]>=vector[lindex])):
            vector[j]=vector[j]+(pert*np.random.randn())*h[j];
    else:
       vector=vectorr[:];
    #print("vector=",vector)
    return  vector
