#In interline.jl
#---Function to interpolate linearly
function interline(uo::Float64,u1::Float64,dy::Int)
Uinterp=zeros(Float64,dy);
const du=u1-uo;
const m=du/dy;
  for i=1:dy
      Uinterp[i]=m*i+uo;
  end

return Uinterp  
end
