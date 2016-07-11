#In halton.jl
using Distributions
function halton(p,n)
Lnd=round(Int,ceil(log(n)/log(p)));                    #--Largest number of digits
b=zeros(Lnd);
u=zeros(n);
uni=zeros(n);
a=-0.5;
c=0.5;
for j=1:n
  i=1;
  b[1]=b[1]+1;
  while b[i]>(p-1+eps())
    b[i]=0;
    i=i+1;
    b[i]=b[i]+1;
  end
  u[j]=0;
  uni[j]=a+(c-a)*u[j];                                  #---Create an uniform distribution
  for k=1:length(b[:])
      u[j]=u[j]+b[k]*p^(float(-k));
      uni[j]=a+(c-a)*u[j];
  end
end
return uni
end  #---end function ran_cauchy
