##In vorticity_sh.jl
function vorticity_sh(u::Array,v::Array,dX::Float64,dY::Float64)
#----Funtion to compute the shear pure contribution of the vorticity field
#----Input parameters
#----u,v Smoothed velocitys
#----dX,dY   lenght differentials of the grid
#----Output  Parameters
#----W vorticity field for the grid m x n
#----Wsh  pure shear contribution of the vorticity field for the grid m x n 
m,n=size(u);
C=ones(m,n)*NaN;
W=ones(m,n)*NaN;
Wsh=ones(m,n)*NaN;                       #---vorticity pure shear contribution
S=ones(m,n)*NaN;                         #---eigenvalues
dudx=0;
dudy=0;
dudyhat=0;
dudysh=0;
dvdx=0;
dvdxhat=0;
dvdxsh=0;
dvdy=0;
kx=0.5/dX;
ky=0.5/dY;
  for i=2:m-1
    for j=2:n-1
      if isnan(u[i,j])==0 && ((isnan(u[i-1,j])==0 && isnan(u[i+1,j])==0)     #--Make sure don't use NaN velocities
                          && (isnan(u[i,j-1])==0 && isnan(u[i,j+1])==0))
         dudx= kx*(u[i+1,j]-u[i-1,j]);                   
         dudy= ky*(u[i,j+1]-u[i,j-1]);                   
         dvdx= kx*(v[i+1,j]-v[i-1,j]);                   
         dvdy= ky*(v[i,j+1]-v[i,j-1]);
         W[i,j]= 0.5*(dvdx-dudy);                                        #--vorticity
         S[i,j]= 0.5*sqrt(4*dudx^2+(dudy+dvdx)^2);
         dudyhat= S[i,j]-W[i,j];                                         #----BRF frame
         dvdxhat= S[i,j]+W[i,j];                                         #----BRF frame
         dudysh= dudyhat-sign(dudyhat)*min(abs(dudyhat), abs(dvdxhat));  #----sh component sh12
         dvdxsh= dvdxhat-sign(dvdxhat)*min(abs(dudyhat), abs(dvdxhat));  #----sh component sh21
         Wsh[i,j] = dvdxsh-dudysh;
      end   
    end   
  end
return W,Wsh
end
