function wrmsfilter(wraw,threshold)
#----Input parameters
#----wraw:         nxm matrix with the raw vorticity field        [s^-1]
#----threshold:    noise level to suppress in the vorticiti field [s^-1]
#----wfiltered:    filtered vorticity field with NaN values under the threshold
m,n=size(wraw);
wfiltered=zeros(m,n);

for i=2:m-1
  for j=2:n-1
     if isnan(wraw[i,j])==0                                    #---- filter just numeric values
        if abs(wraw[i,j])<(threshold)                          #---- threshold<wraw<threshold
           wraw[i,j]=NaN;
        end
     end   
  end
end
  
wfiltered=wraw;
return wfiltered
end
