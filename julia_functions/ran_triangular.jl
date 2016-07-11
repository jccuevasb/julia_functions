#In ran_uniform.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
function ran_triangular(a,b,c)

  U= rand();
  f= (c-a)/(b-a);
  if U<=f
     rnum= a+sqrt(U*(b-a)*(c-a));
  else
     rnum= b-sqrt((1-U)*(b-a)*(b-c));
  end   

return rnum
end  #---end function ran_cauchy
