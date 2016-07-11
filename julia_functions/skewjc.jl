#In skewjc.jl
include("meanjc.jl")                          #--include library 
include("varjc.jl")                          #--include library 

function skewjc(vector::Array)
N=length(vector);
s3=0.0;                                        #--sample mean sum
m=0.0;                                        #--counter for the sample size
mu=meanjc(vector);                           #---mean

#print("mu= ",mu);
#print("\n");
sigma2=varjc(vector);                        #---variance
#print("sigma2= ",sigma2);
 if sigma2==0;
    skewness=0;
 elseif isnan(sigma2)==1;
    skewness=NaN;
 else   
   for i=1:N
      if isnan(vector[i])==0;                  #--make sure to sum numeric values
         s3=s3+vector[i]^3;
         m=m+1;
      end
   end
   skewness=(s3/m-3*mu*sigma2-mu^3)/sigma2^(3/2);
 end  
return skewness
end
