#In ran_gamma.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Distributions
#import numpy as np
function ran_gamma(shape,scale)
#scale parameter beta=1/theta, shape parameter alpha=k
   if shape <= 0.0
      error("Shape parameter must be positive")
   end
   if scale <= 0.0
      error("Scale parameter must be positive")
   end
                          
 ## Implementation based on "A Simple Method for Generating Gamma Variables"
 ## by George Marsaglia and Wai Wan Tsang.  
 ## ACM Transactions on Mathematical Software
 ## Vol 26, No 3, September 2000, pages 363-372.
   if shape >= 1.0
      d = shape - 1.0/3.0
      c = 1.0/sqrt(9.0*d)
      while true
         x = randn();                     #---rand normal 0 mean 1 std
         v = 1.0 + c*x
         while v <= 0.0
             x = randn()
             v = 1.0 + c*x
         end
         v = v*v*v
         u = rand()
         xsq = x*x
         if u < 1.0 -.0331*xsq*xsq || log(u) < 0.5*xsq + d*(1.0 - v + log(v))
            return scale*d*v
         end
      end
   else
      g = ran_gamma(shape+1.0, 1.0)
      w = rand()
      #return scale*g*pow(w, 1.0/shape)
      return scale*g*w^(1.0/shape)
  end    
end  #---end function ran_gamma
