#In meanprod.jl
function meanprod(vector1,vector2,ydata)
N=length(vector1);
svv=0;                                               #--sample mean sum
m=0;                                                 #--initial size of the sample mean
 for i=1:N
   if (isnan(vector1[i])==0 && isnan(vector2[i])==0) && (isnan(ydata[i])==0)
       m=m+1;                                        #---counter for the sample size
       svv=svv+(vector1[i]*vector2[i]);
   end
 end
mu=svv/m;
return mu
end
