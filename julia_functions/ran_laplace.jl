#In ran_laplace.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
function ran_laplace(mean,scale)
  if scale <= 0.0
     error("Scale parameter must be positive")
  end
u = rand();
  if u < 0.5
    retval = mean + scale*log(2.0*u) 
#    print("u= ",u)
  else
    retval = mean - scale*log(2*(1-u))
  end
return retval
end  #---end function ran_cauchy
