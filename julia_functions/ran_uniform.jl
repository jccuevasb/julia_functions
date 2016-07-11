#In ran_uniform.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
function ran_uniform(a,b)
   rnum=a+(b-a)*rand();
return rnum
end  #---end function ran_cauchy
