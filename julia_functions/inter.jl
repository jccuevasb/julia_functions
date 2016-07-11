#in inter.jl
using Interpolations
include("interp1jc.jl");

function   inter(y1::Array{Float64,1},y2::Array{Float64,1},yplus::Array{Float64,1})

cp=diff(sign(y1-y2));
xpos=0;
ypos=0;
x=Float64[];
  if sum(cp)!=0                     #--summ=0 no intersection points
     indey1=find(g-> g!=0,cp)
     if indey1[length(indey1)] < length(yplus)
        p1=y1[indey1];
        p2=y1[indey1+1];
        p=[p1' ; p2'];                                       #----create a matrix
        x1=yplus[indey1];
        x2=yplus[indey1+1];
        x=[x1' ; x2'];                                       #----create a matrix
     else
        p1=y1[indey1];
        p2=0;
        p=[p1' ; p2'];
        x1=yplus[indey1];
        x2=0;
        x=[x1' ; x2'];
     end

     interval=(x1+x2)/2;

     if ndims(x)==1      #--check if it is an array or vector
        m=length(x);     #--rows
        n=1;             #--column
     else
        m,n=size(x);
     end   

     px=zeros(n);

     for i=1:n           #--iterate over first dimmension of tuple m=(row,column)
         itp=interpolate((vec(x[:,i]),),vec(p[:,i]),Gridded(Linear()));
         px[i]=itp[interval[i]];
     end

     C,index2=findmax(abs(px));
     xpos=interval[index2];
     ypos=px[index2];

  else
  xpos=0;
  ypos=0;
  end   

return xpos,ypos
end
