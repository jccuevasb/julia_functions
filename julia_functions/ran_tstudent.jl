#In ran_uniform.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
function ran_tstudent(dof)
 
  if dof <= 0
     error("Degrees of freedom must be positive")
  end

  ## See Seminumerical Algorithms by Knuth
  y1 = randn();
  y2 = rand_chi_square(dof);
  rnum = y1 / sqrt(y2 / dof);

return rnum
end  #---end function ran_cauchy
