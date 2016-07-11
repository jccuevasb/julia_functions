#In meanjc.jl

function meanjc(vector)
N=length(vector);
s=0.0;                                        #--sample mean sum
m=0.0;                                        #--counter for the sample size
  for i=1:N
     if isnan(vector[i])==0;                  #--make sure to sum numeric values
        m=m+1;
        s=s+vector[i];
     end
  end
mu=s/m;
mu
end
