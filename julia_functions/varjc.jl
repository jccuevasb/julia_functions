#In varjc.jl
#https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
include("meanjc.jl")                          #--include library 

function varjc(vector)
N=length(vector);
ss=0.0;                                        #--sample square mean sum
s=0.0;                                        #---sample mean sum
m=0.0;                                        #--counter for the sample size
  for i=1:N
     if isnan(vector[i])==0;                  #--make sure to sum numeric values
        K=vector[1];                          #--variance is invariant respect to changes in a location parameter
        s=s+vector[i]-K
        ss=ss+(vector[i]-K)*(vector[i]-K);
        m=m+1;
     end
  end
#sigma2=(ss/m)-meanjc(vector)^2;                #---This approach gives rounding errors and negative values
sigma2=(ss-(s^2/m))/m;                #---
# use n instead of (m-1) if want to compute the exact variance of the given data
# use (m-1) if data are samples of a larger population
return sigma2
end
