#In interp1jc.jl
#Function to interpolate data sets linearly
function interp1jc(x,y,grid)
n=length(grid);
yp=zeros(n);
m=0.0;                             #--Slope
for i=1:length(x)-1
m=(y[i+1]-y[i])/(x[i+1]-x[i]);
  for j=1:n
    if grid[j]>x[i] && grid[j]<x[i+1]
       yp[j]=y[i]+(grid[j]-x[i])*m;
    elseif grid[j]==x[i]
       yp[j]=y[i];
    elseif grid[j]==x[i+1]
       yp[j]=y[i+1];
    end
  end
end
return yp
end
