#In nankurt.jl
#https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
function nankurt(vector::Array)
N=length(vector);  
m=0.0;             #---Initialize counter for samples
mu=0.0;            #---Initialize mean
M2=0.0;            #---Intitialize sum(x_i-mu)^2
M3=0.0;            #---Intitialize sum(x_i-mu)^2
M4=0.0;            #---Intitialize sum(x_i-mu)^2

 for i=1:N
      if isnan(vector[i])==0;                  #--make sure to sum numeric values
         n1=m;
         m=m+1;
         delta=vector[i]-mu;
         delta_m=delta/m;
         delta_m2=delta_m*delta_m;
         term1=delta*delta_m*n1;
         mu=mu+delta_m;
         M4=M4+term1*delta_m2*(m*m-3*m+3)+6*delta_m2*M2-4*delta_m*M3
         M3=M3+term1*delta_m*(m-2)-3*delta_m*M2
         M2=M2+term1;
      end
 end
#kurtosis=(m*M4)/(M2*M2)-3;      #--normalized
kurtosis=(m*M4)/(M2*M2); 
return kurtosis
end
