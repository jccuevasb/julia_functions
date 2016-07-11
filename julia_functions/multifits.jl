# In multifits.jl
#--This function computes multiples fits for dat
using Optim       #load optimization package
function multifits(xdata::Array,ydata::Array,strfit::ASCIIString)
#--Input parameters
#--x= data values in x coordinate (row-vector)
#--y= data values in y coordinate (row-vector)
#--strfit= type of fit 'logarithmic' or 'linear'
#--Output parameters
#--coeff=coefficients for the fit i.e. logarithmic, linear, etc
#--fminres=value of the fit in the given coefficients

function expfunc(betas)  #---exponetial fit
         pred_i = betas[1].*exp(-betas[2]*xdata);
         err=sum((ydata - pred_i).^2);
     return err
end

function logfunc(betas)  #---logarithmic fit
         pred_i = betas[1]*log(xdata)+betas[2];    #y=A*log(x)+B
         err=sum((ydata - pred_i).^2);
     return err
end

function linearfunc(betas)  #---linear fit
         pred_i = betas[1]*xdata+betas[2];         #y=A*x+B
         err=sum((ydata - pred_i).^2);
     return err
end

function powfunc(betas)  #---power fit
         pred_i = betas[1]*xdata.^(betas[2]);      #y=A*x^B
         err=sum((ydata - pred_i).^2);
     return err
end

start=[0.0,0.5];   #initial guess
if ((strfit=="exponential")==1)
   res = optimize(expfunc,start, method = :nelder_mead);
elseif  ((strfit=="logarithmic")==1)
   #res = optimize(logfunc,start, method = :nelder_mead);
   res = optimize(logfunc,start, method = Optim.NelderMead());
elseif  ((strfit=="linear")==1)
   res = optimize(linearfunc,start, method = :nelder_mead);
elseif  ((strfit=="power")==1)
   res = optimize(powfunc,start, method = :nelder_mead);
end

coeff=res.minimum;          #coefficients of minimum
fminres=res.f_minimum;      # Value of the function at minimum
return coeff, fminres
end
