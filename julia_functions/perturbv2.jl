#In perturbv2.jl

function perturbv2(yplinear::Vector,Np::Integer)
inistep=Int64[];
lstep=Int64[];
pert=Float64[];
peru=Float64[];
Nve=length(yplinear);
ypert=zeros(Nve);
vector=zeros(Nve);      #---Initialize vector to storage perturbed positions
vector=log(yplinear);
pert=1.0;
peru=1.3;
#print("ypert=",ypert)
 if Np>0;
    inistep=2;                          #--initial step to perturb;
    lstep=Nve-1;                        #--Last step to perturb
    h=zeros(lstep);                     #--Initialize height vector
    uniper=zeros(lstep);                #--Initialize uniform perturbation
#    print(" vector= ",vector);
    for j=inistep:lstep                       #--iterate from the 2nd yp until the Nth-1 yp step
      #uniper[j]=-peru+(2*peru)*rand();         #---compute percentage of uniform perturbation
      uniper[j]=pert*randn();                  #---compute percentage of gaussian perturbation
      #print("uniper= ",uniper[j])
      if uniper[j]<0;                          #>0 positive biased    <0 negative biased
          h[j]=vector[j+1]-vector[j];         #--Actual height of the y(j) step  y[j+1]-y[j]  upper step
      else
          h[j]=vector[j]-vector[j-1];            #--Actual height of the y(j) step  y[j]-y[j-1]  lower step
      end
      ypert[j]=vector[j]+uniper[j]*h[j];                           #--uniform  (-1,1)
      while ((ypert[j]<=vector[1]) | (ypert[j]>=vector[Nve]))
           #uniper[j]=-peru+(2*peru)*rand();         #---compute percentage of uniform perturbation
           uniper[j]=pert*randn();                  #---compute percentage of gaussian perturbation
           if uniper[j]<0;
              h[j]=vector[j+1]-vector[j];         #--Actual height of the y(j) step  y[j+1]-y[j]  upper step
           else
              h[j]=vector[j]-vector[j-1];            #--Actual height of the y(j) step  y[j]-y[j-1]  lower step
           end  
           #vector[j]=vector[j]+uniper[j]*h[j];                           #--uniform  (-1,1)
           ypert[j]=vector[j]+uniper[j]*h[j];                           #--uniform  (-1,1)
      end            
    end         #----end for loop
    ypert[1]=vector[1];                                    #---first yp position
    ypert[Nve]=vector[Nve];                                #---last yp position
    #print("ypert=",ypert)
    #print("\n");
 else   
    ypert=vector;
 end             #----end if-else

ypert=exp(ypert);

return ypert
end
