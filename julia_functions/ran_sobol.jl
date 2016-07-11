#In ran_sobol.jl
#---This is a function to generate random numbers that follow a cauchy distribution
#---Taken from http://www.johndcook.com/julia_rng.html
using Sobol
#function ran_sobol(Np::Integer)
function ran_sobol(s)
#---use the enthropy of the system to generate random uniform distributed numbers
#dim=1;
#rnum=0;
#s=SobolSeq(dim);                            #---A sequence of numbers
test=next(s);                               #---An array of one value
#print("seq= ",s);
rnum=test[1];
#rnum=Np+5;
return rnum
end  #---end function ran_cauchy
