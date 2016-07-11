#In com_peaks.jl
function com_peaks(phi_uu::Array{Float64,2},yplus::Array)
m,n=size(phi_uu);                  #Number of columns
#----Find max value per columns 
#C,index=max(abs(phi_uu),[],1);    #--C values, I index value in Phi_uu
C,index=findmax(abs(phi_uu),1);    #--C values, I index value in Phi_uu
#C=C[:];                            #--convert to row vector
#index=index[:];                    #--convert to row vector
#indexcol=[0:m:n*m-m];
indexcol=collect(0:m:n*m-m);
index=index'-indexcol;
max_val=zeros(n);
ycoord=zeros(n);
 for i=1:n
     max_val[i]=phi_uu[index[i],i]; #-Create an array with the max values
     ycoord[i]=yplus[index[i]];     #-Create an array with the yplus pos
 end
return max_val,ycoord
end
