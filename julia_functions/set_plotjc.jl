#In set_plotjc.jl
#---This script attempts to create a really nice plots
#---for publication
#Using PyPlot
#using PyCall
function set_plotjc(ax)
fs=17;
w=1.4;         #---width of the ticks
setp(ax[:get_yticklabels](),color="black",fontsize=fs) # Y Axis font formatting
setp(ax[:get_xticklabels](),color="black",fontsize=fs) # Y Axis font formatting
#########################
#  Set tick dimensions  #
#########################
ax[:xaxis][:set_tick_params](which="major",length=8,width=w)    #--length=hight   xaxis
ax[:yaxis][:set_tick_params](which="major",length=8,width=w)    #--length=hight   yaxis
ax[:xaxis][:set_tick_params](which="minor",length=5,width=w)
ax[:yaxis][:set_tick_params](which="minor",length=5,width=w)
#PyDict(pyimport("matplotlib")["rcParams"])["font.sans-serif"] = ["Courier"]
return ax
end
