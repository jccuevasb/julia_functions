#In nanvar.jl
#https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Higher-order_statistics
function nanvar(vector::Array)
N=length(vector);  
m=0.0;             #---Initialize counter for samples
mu=0.0;            #---Initialize mean
M2=0.0;            #---Intitialize sum(x_i-mu)^2

 for i=1:N
      if isnan(vector[i])==0;                  #--make sure to sum numeric values
         m=m+1;
         delta=vector[i]-mu;
         mu=mu+delta/m;
         M2=M2+delta*(vector[i]-mu);
      end
 end
 if m<2;
    return NaN
 else
    return M2/m
 end   
end
