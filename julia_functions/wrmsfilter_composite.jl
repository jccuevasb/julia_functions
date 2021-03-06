function wrmsfilter_composite(wsh::Array{Float64,2},y::Array{Float64,2},wshmean::Array{Float64,2},wshrms::Array{Float64,2})
#----Input parameters
#----wsh:          nxm matrix with the raw vorticity sh field        [s^-1]
#----wshrms:       noise level to suppress in the vorticiti field [s^-1]
#----wshmean:      mean of pure shear vorticity contribution
#----y:            wall normal position matrix
#----Output parameters
#----wfiltered:    filtered vorticity field with NaN values under the threshold
m,n=size(wsh);
wfiltered=zeros(m,n);
ywall=5.0;       #[mm]
cut=1;
for i=2:m-1
  for j=2:n-1
     if isnan(wsh[i,j])==0                                    #---- filter just numeric values
        if y[i,j]<=ywall                                       #---- composite filter in the wall
           if (abs(wsh[i,j])<=abs(wshmean[i,j])) | (abs(wsh[i,j])<=cut*wshrms[i,j])
              wsh[i,j]=NaN;
           end   
        #elseif abs(wsh[i,j])<=cut*wshrms[i,j]                #----free-stream filter
        elseif abs(wsh[i,j])<=210                #----free-stream filter
            wsh[i,j]=NaN;
        end    
     end   
  end
end
  
wfiltered=wsh;
return wfiltered
end
