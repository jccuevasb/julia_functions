#in lsqf_smooth.jl
#script created by follow the steps in
#http://lls.readthedocs.org/en/latest/julia_examples.html
function lsqf_smooth(ynoise,lambda)
#This function satisfy the following objective
#sum (s_i-x_i)^2+lambda*sum (s_i-s_i-1)^2         i=1,...n
#input parameters
#ynoise= vector to be smoothed
#lambda= smooth parameter  Larger lambda, smoother model will be.
n=length(ynoise);
ydata=zeros(n);
smooth_objective = sum((ydata[1:n-1]-ydata[2:n]).^2);
#ysmooth=Variable(n);            # A vector variable with n rows and 1 column




return n
end

