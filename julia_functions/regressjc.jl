#In regressjc.jl
include("meanprod.jl")
include("varjc.jl")

function regressjc(y,X)
#---function to compute the linear regression of multiple models
#---inputs; y= true data (dependent variable )time series column vector
#---x=matrix formed by column vectors of predictor variables
#---X=(x1;x2,x3,...,xn)
#---Outputs:
#---alphas= regression coefficients, yhat=model estimate, shat= skill estimate.
#---this function satisfies the model
#---yi=alpha_1*x_i1+alpa_2*x_2i+....+alpha_n*x_in
m,n=size(X);                                     #--compute size of the matrix
#--Initialize matrix X and vector solution zm Dmm*alphas_m=zm
Dmm=zeros(n,n);
zm=zeros(n);                                     #--<yxm>
shat=0;
 for i=1:n
   for j=1:n
     Dmm[i,j]=meanprod(X[:,i],X[:,j],y);
   end
   zm[i]=meanprod(y,X[:,i],y);
 end
alphas=Dmm\zm;                                   #--alphas_m=Dmm^-1*zmy
#---compute the estimate model y_hat=sum aplha_m*xm
yhat=zeros(m);
 for k=1:n;
   yhat=yhat+(alphas[k]*X[:,k]);
 end
shat=varjc(yhat)/varjc(y);                       #--skill estimate model sigma^2_yhat/sigma^2_ydata
return alphas, yhat, shat
end
