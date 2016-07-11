##In vorticity_8points.jl
function vorticity_8points(u::Array,v::Array,dX::Float64,dY::Float64)
#--function that return the vorticity based on the eight-points approximation
#--see [1] Adrian & Whesterweel, Particle Image Velocimetry, Cambridge University Press, 2011
#--Compute the circulation around the point i,j using the information
#--of the eight neighboring points
m,n=size(u);
C=ones(m,n)*NaN;
W=ones(m,n)*NaN;
  for i=2:m-1
    for j=2:n-1
      if isnan(u[i,j])==0 && ((isnan(u[i-1,j+1])==0 && isnan(u[i+1,j+1])==0)
                          && (isnan(u[i-1,j-1])==0 && isnan(u[i+1,j-1])==0))
       C[i,j] = 0.5*dX*(u[i-1,j-1] + 2*u[i,j-1] + u[i+1,j-1])+
                0.5*dY*(v[i+1,j-1] + 2*v[i+1,j] + v[i+1,j+1])-
                0.5*dX*(u[i+1,j+1] + 2*u[i,j+1] + u[i-1,j+1])-
                0.5*dY*(v[i-1,j+1] + 2*v[i-1,j] + v[i-1,j-1]);
       # vorticity
       W[i,j] = C[i,j]/(4*dX*dY);
       end   
    end   
  end
return W  
end
