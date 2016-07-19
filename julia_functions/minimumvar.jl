#---In minimumvar.jl
#---Function to compute the minimu value of the wzl/wztotal ratio
#---and its position
function minimumvar(ycoord,ratio)
va,ind=findmin(ratio);
ypl=ycoord[ind];

return ypl,va,ind
end
