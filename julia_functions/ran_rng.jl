#In ran_uniform.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
#using Distributions
function ran_rng(Np::Integer)
#---use the enthropy of the system to generate random uniform distributed numbers
#rng = MersenneTwister();      #  Initialize the default RNG!!!  
#rng=  RandomDevice();          #--It is limited if you have many iterations
#srand(Np);
srand();                       #--It is limited also if you have many iterations ()use OS clock to produce a seed
rnum=rand();
return rnum
end  #---end function ran_cauchy
