#In perturb.jl
#include("ran_beta.jl")
#include("ran_triangular.jl")
#include("ran_rng.jl")
#include("ran_laplace.jl")
#using PyCall

function perturb(vector::Vector,Np::Integer,wvf::Integer)
#----wvf width of vortical fissure
#@pyimport pylab as npy
inistep=Int64[];
lstep=Int64[];
Nve=length(vector);
ypert=zeros(Nve);      #---Initialize vector to storage perturbed positions
const pert=1.0;
const peru=1.5;
alpha=2;           #----parameters for beta distribution
beta=2;            #----Parameters for beta distribution
 if Np>0;  
    inistep=2;                          #--initial step to perturb;
    lstep=Nve-1;                        #--Last step to perturb Nth-1
    #lstep=Nve;                           #--Last step to perturb Nth
    h=zeros(lstep);           #--Initialize height vector
#    print(" vector= ",vector);
    for j=inistep:lstep                       #--iterate from the 2nd yp until the Nth-1 yp step
      h[j]=vector[j]-vector[j-1];            #--Actual height of the y(j) step
      ypert[j]=vector[j]+(pert*randn()*h[j]);     #--gauss distr mu=0 sigma=1 (-5,5) normalized by 5 (check this for problems in b.c)
      #ypert[j]=vector[j]+(ran_laplace(0,0.5)*h[j]);     #--laplace distr mu=0 scale factor=0.5
      #ypert[j]=vector[j]+(pert*randn()*h[j])+40*randn();         #--gauss distr mu=0 sigma=1 
      #vector[j]=vector[j]+(pert*randn())*h[j]+(-peru+(2*peru)*rand());              #--gauss distr + uniform
      #ypert[j]=vector[j]+(-peru+(2*peru)*rand())*h[j];                           #--uniform  (-1,1)
      #ypert[j]=vector[j]+ran_triangular(-0.5,0.5,0)*h[j];                           #--triangular  (-0.5,0.5)
      #ypert[j]=vector[j]+(-peru+(2*peru)*ran_rng(Np))*h[j];                           #--use rng device to produce random numbers
      #ypert[j]=vector[j]+ran_beta(alpha,beta)*h[j];                                       #--beta  (-0.5,0.5)
      #ypert[j]=vector[j]+(-peru+(2*peru)*rand())*h[j]+(pert*randn());             #--uniform + gaussian
      #if j==2;                                                           #---constraint second vf
      #    while ((ypert[j]<=vector[1]) | (ypert[j]>=vector[3]))
          #while ((ypert[j]<=vector[1]) | (ypert[j]>=vector[3])) && ((ypert[j]-ypert[j-1])<7)
      #         ypert[j]=vector[j]+(pert*randn()*h[j]);   
      #    end 
      #end
      while ((ypert[j]<=vector[1]) | (ypert[j]>=vector[Nve]))  #--Constraint the first and last step
      #while ((ypert[j]<=(vector[1]+wvf)) | (ypert[j]>=(vector[Nve]-wvf)))  #--Constraint the first and last step
      #while ((ypert[j]<=vector[1]) | (ypert[j]>=vector[Nve])) && ((ypert[j]-ypert[j-1])<7)  #--Constraint the first and last step
      #while (ypert[j]<vector[1])                                                  #--Constraint first step <0
           ypert[j]=vector[j]+(pert*randn()*h[j]);                    #--gauss +  gauss
           #ypert[j]=vector[j]+(ran_laplace(0,0.5)*h[j]);     #--laplace distr mu=0 scale factor=0.5
           #ypert[j]=vector[j]+(pert*randn()*h[j])+40*randn();         #--gauss distr mu=0 sigma=1 
           #ypert[j]=vector[j]+(pert*randn())*h[j]+(-peru+(2*peru)*rand());         #--gauss distr + uniform
           #ypert[j]=vector[j]+(-peru+(2*peru)*rand())*h[j];                           #--uniform  (-1,1)
           #ypert[j]=vector[j]+ran_triangular(-0.5,0.5,0)*h[j];                           #--triangular  (-0.5,0.5)
           #ypert[j]=vector[j]+(-peru+(2*peru)*ran_rng(Np))*h[j];                           #--use rng device to produce random numbers
           #ypert[j]=vector[j]+ran_beta(alpha,beta)*h[j];                                       #--beta  (-0.5,0.5)
           #print("ypert[j]=",ypert[j])
           #print("\n");
           #ypert[j]=vector[j]+(-peru+(2*peru)*rand())*h[j]+(pert*randn());             #--uniform + gaussian
      end            
    end
    ypert[1]=vector[1];                                    #---first yp position
    ypert[Nve]=vector[Nve];                                #---last yp position
 else
    ypert=vector;
 end
       

return ypert
end
