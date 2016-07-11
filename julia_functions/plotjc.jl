#In plotjc.jl
#---This script attempts to create a really nice plots
#---for publication
#Using PyPlot
function plotjc(xdata::Array, ydata::Array, col, marker, lw)
fs=17;
w=1.4;         #---width of the ticks
h1=plot(xdata,ydata,color=col,markersize=3,marker,markeredgecolor=col,markerfacecolor="none",linewidth=lw);
ax=gca();
setp(ax[:get_yticklabels](),color="black",fontsize=fs) # Y Axis font formatting
setp(ax[:get_xticklabels](),color="black",fontsize=fs) # Y Axis font formatting
#########################
#  Set tick dimensions  #
#########################
ax[:xaxis][:set_tick_params](which="major",length=8,width=w)    #--length=hight   xaxis
ax[:yaxis][:set_tick_params](which="major",length=8,width=w)    #--length=hight   yaxis
ax[:xaxis][:set_tick_params](which="minor",length=5,width=w)
ax[:yaxis][:set_tick_params](which="minor",length=5,width=w)
#setp(ax[:get_ylabel](),labelsize=fs) # Y Axis font formatting
#setp(ax[:get_ylabel](),fontsize=fs) # Y Axis font formatting
return h1
end
