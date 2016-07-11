#In ran_uniform.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
function rand_chi_square(dof)

  rnum = ran_gamma(0.5*dof, 2.0);

return rnum
end  #---end function ran_cauchy
