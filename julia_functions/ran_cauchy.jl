#In ran_cauchy.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
function ran_cauchy(median,scale)
   if scale <= 0.0
      error("Scale parameter must be positive")
   end
rn=rand(Cauchy(median,scale));
return rn

end  #---end function ran_cauchy
