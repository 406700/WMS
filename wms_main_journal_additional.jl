
include("wms_main_functions.jl")

# using PyCall
# pyimport("matplotlib")
# plt = pyimport("matplotlib.pyplot")
# mpl = pyimport("matplotlib")
# mpl.use("Agg")  # Set the backend to Agg

# Now import PyPlot
# using PyPlot

rkey=Dict(zip(["λ", "ν", "L", "L_prime", "L_direct_avg", "L_prime_direct","FP_array", "LD_array","scale_factor","Δim"],[1:10;]))

results=load_results("20240110__combined_results")
results2=load_results("20240110__combined_results_17")

λ_17=1546.546e-9
λ_26=1546.637e-9

modulation_values=[5,10,15,20,25,30]

scaling_factors=[results[value][rkey["scale_factor"]] for value in modulation_values]
scaling_factors2=[results2[value][rkey["scale_factor"]] for value in modulation_values]

#confirm the scaling Plots.jl
# p=plot()
# for (index,value) in enumerate(modulation_values)
#     plot!(p,results[value][rkey["L_prime"]].*scaling_factors[index])
#     plot!(p,results[value][rkey["L_prime_direct"]],label="direct")
# end
# display(p)

######################################################################################scale factors

Δν_h=[1/x for x in scaling_factors]*-3e10 #L_0'=scale_factor/3e-10 *L'    -3e10/scale_factor= Δν<h>
Δν_h2=[1/x for x in scaling_factors2]*-3e10 #
c=2.99e8


fig, ax = plt.subplots()
ax.set_ylabel(L"$|\Delta \lambda \langle h \rangle |$")
ax.set_xlabel("modulation voltage (mV)")  # Optionally enforce fontsize
ax.plot(modulation_values,abs.(Δν_h*λ_26^2/c*1e12),label="26 mA",marker=".")
ax.plot(modulation_values,abs.(Δν_h2*λ_26^2/c*1e12),label="17 mA",linestyle="--",marker=".")
plt.legend(loc="upper left")
# Show the plot
plt.tight_layout(pad=0.4)
plt.savefig("./matplotlibfigures/scaling_factors.png")

# p=plot(xlabel="modulation voltage (mV)",ylabel="|Δλ<h>|",legend=:bottomright,xlims=(0, maximum(modulation_values)),ylims=(0,35),linestyle=:solid, linealpha=0.5, linewidth=4*upscale,left_margin=5mm*upscale,grid=false)
# plot!(p,modulation_values,abs.(Δν_h*λ_26^2/c*1e12),label="26 mA")
# plot!(p,modulation_values,abs.(Δν_h2*λ_26^2/c*1e12),label="17 mA")
# #savefig("./figures/scaling_factors.svg")
# p=empty 

#functions for centering
##############################################################################load the data and center
centered_results=center_pairs_in_dict(results,modulation_values)
centered_results=center_lambda(centered_results,modulation_values)
results=empty

#centered_results2=center_pairs_in_dict(results2,modulation_values)
#centered_results2=center_lambda(centered_results2,modulation_values)
results2=empty 
GC.gc() #clear memory

##########################################################################plot the scaled data for 26mA and add inset.
#


fig, ax = plt.subplots()
ax.set_ylabel("T")
ax.set_xlabel("wavelength shift (pm)")  # Optionally enforce fontsize
for (index,value) in enumerate(modulation_values)
    flab=modulation_values[index]
    ax.plot(centered_results[value][1]*1e12,centered_results[value][rkey["L"]],label="$flab mV")
end
ax.plot(centered_results[5][1]*1e12,centered_results[5][rkey["L_direct_avg"]],label=L"$T_0$")
plt.legend(loc="upper left")
# Show the plot
plt.tight_layout(pad=0.4)
plt.savefig("./matplotlibfigures/L.png",dpi=800)

fig, ax = plt.subplots()
ax.set_ylabel(L"$T'$ (1/Thz)")
ax.set_xlabel("modulation voltage (mV)")  # Optionally enforce fontsize
for (index,value) in enumerate(modulation_values)
    flab=modulation_values[index]
    ax.plot(centered_results[value][1]*1e12,centered_results[value][rkey["L_prime"]].*scaling_factors[index]*1e12,label="$flab mV")
end
ax.plot(centered_results[5][1]*1e12,centered_results[5][rkey["L_prime_direct"]]*1e12,label=L"$T_0'$")#"\$T\\prime_0\$"
plt.legend(loc="upper left")
# Show the plot
plt.tight_layout(pad=0.4)
plt.savefig("./matplotlibfigures/L_prime.png",dpi=800)

#p=plot(xlabel="wavelength shift (pm)",ylabel=L"L",linestyle=:solid, linealpha=0.5, linewidth=4*upscale,frame=:box,left_margin=5mm*upscale,grid=false)
# for (index,value) in enumerate(modulation_values)
#     flab=modulation_values[index]
#     plot!(p,centered_results[value][1]*1e12,centered_results[value][rkey["L"]],label="$flab mV")
#     #plot!(p,results[value][rkey["L_prime_direct"]],label="direct")
# end
# plot!(p, centered_results[5][1]*1e12,centered_results[5][rkey["L_direct_avg"]],label=L"L_0")
# display(p)
#savefig("./figures/L")

# p2=plot(xlabel="wavelength shift (pm)",ylabel="T' (1/THz)",linestyle=:solid, linealpha=0.5, linewidth=4*upscale,frame=:box,left_margin=5mm*upscale,grid=false)
# for (index,value) in enumerate(modulation_values)
#     flab=modulation_values[index]
#     plot!(p2,centered_results[value][1]*1e12,centered_results[value][rkey["L_prime"]].*scaling_factors[index]*1e12,label="$flab mV")
# end
# plot!(p2, centered_results[5][1]*1e12,centered_results[5][rkey["L_prime_direct"]]*1e12,label="T₀'")
# display(p2)
#savefig("./figures/L_prime")

################################################################################## add an inset plot 
# # Define the range for the inset plot [t1:t2]
# t1 = 350# Define t1
# t2 = 450 # Define t2

# # Create the inset plot
# inset_p2 = plot(legend=false)  
# for (index, value) in enumerate(modulation_values)
#     plot!(inset_p2, centered_results[value][rkey["L_prime"]][t1:t2].*scaling_factors[index])
#     # Add other series if needed, make sure to restrict to [t1:t2]
# end
# display(inset_p2)
# # Position the inset plot within the main plot
# # The arguments are (plot, [left, bottom, width, height]), values range from 0 to 1
# plot!(centered_results[5][rkey["L_prime"]][t1:t2].*scaling_factors[1],
#     frame=:box,
#     yticks=false,
#     xticks=false,
#     inset=bbox(0.5, 0.5, 0.3, 0.3),
#     subplot=2
# )
# display(p2)
#check the error compared to the direct derivative

################################################################################plot the derivative error
# scaling=1/0.8
# fntsm = Plots.font("sans-serif", pointsize=round(12.0*upscale)*scaling)
# fntlg = Plots.font("sans-serif", pointsize=round(16.0*upscale)*scaling)


modulation_values=[10,20,30]
scaling=scaling_factors[[2,4,6]]
div_error=Dict()

fig, ax = plt.subplots()
ax.set_ylabel(L"$L^\prime-L_0^\prime$ (1/THz)")
ax.set_xlabel("wavelength shift (pm)")  # Optionally enforce fontsize
for (index,value) in enumerate(modulation_values)
    flab=modulation_values[index]
    er=centered_results[value][4]*scaling[index].-centered_results[value][6]
    div_error[value]=er
    ax.plot(centered_results[5][1]*1e12,div_error[value]*1e12, label="$flab mV")
end
plt.legend(loc="upper left")
# Show the plot
plt.tight_layout(pad=0.4)
plt.savefig("./matplotlibfigures/derivative_errors",dpi=800)


# modulation_values=[10,20,30]
# scaling=scaling_factors[[2,4,6]]

# p=plot(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm, xlabel="wavelength shift (pm)",ylabel=L"L^\prime-L_0^\prime"* " (1/THz)",linestyle=:solid, linealpha=0.5, linewidth=4*upscale,frame=:box,left_margin=5mm*upscale,grid=false)
# div_error=Dict()
# for (index, value) in enumerate(modulation_values)
#     flab=modulation_values[index]
#     er=centered_results[value][4]*scaling[index].-centered_results[value][6]
#     div_error[value]=er
#     plot!(p,centered_results[5][1]*1e12,div_error[value]*1e12, label="$flab mV")
# end
# display(p)
#savefig("./figures/derivative_errors.png")

################################################################################plot the demonstration data 
#copy of code in main functions, to make sure it isn't overridden somehow. 
# gr()
# dpi=600
# default(dpi = 600)
# # Width in inches
# width_in_inches = 3.5
# width_in_pixels = width_in_inches * dpi

# # Set width and allow height to scale
# default(size = (width_in_pixels, width_in_pixels))
# upscale = 1 #8x upscaling in resolution
# fntsm = Plots.font("Helvetica", pointsize=round(12.0*upscale))
# fntlg = Plots.font("Helvetica", pointsize=round(16.0*upscale))
# default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
# default(size=(width_in_pixels,width_in_pixels)) #Plot canvas size

# centered_results=empty
# GC.gc()
# #############get data
string="data_config_"*"20"*"ma.jl"
include(string)
FP,time=load_and_return_raw_data(fp_file,ld_file,fp_direct_file,ld_direct_file)
####################plots

gr()

inset_time=0.55

p=plot(xlabel="time (S)",ylabel="Signal (V)",linestyle=:solid, linealpha=0.5, linewidth=4*upscale,frame=:box,legend=false,left_margin=5mm,bottom_margin=10mm*upscale,grid=false)
plot!(p,time[1:100:end],FP[1:100:end])
# #savefig(p,"./figures/raw_data_1.png")

vline!([inset_time],label="")
p2=plot(xlabel="time (S)",ylabel="Signal (V)",linestyle=:solid, linealpha=0.5, linewidth=4*upscale,frame=:box,legend=false,left_margin=10mm,bottom_margin=10mm*upscale,grid=false)
idx=findfirst(x->x>inset_time,time)
plot!(p2,time[idx:1:idx+200000],FP[idx:1:idx+200000],label="")
# #savefig(p2,"./figures/raw_data_2.png")


combined_plot = plot(p, p2, layout=(2, 1))
#savefig("./figures/raw_data.ps")


l = @layout [a b]
ylabel!(p2,"")
combined_plot = plot(p, p2, layout=2,legend=false,link=:y,size=(600*3.5,(600*3.5)/2))
# savefig("./figures/raw_data_sidebyside.png")

#######################################################3#using matplotlib from PyPlot
#
# Set the DPI (resolution)
#overide default aspect ratios
total_width_in_inches = 3.5
num_subplots = 2
subplot_width_in_inches = total_width_in_inches / num_subplots
subplot_height_in_inches = subplot_width_in_inches#aspect_ratio

fig, axs = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))

# Load the data (assuming the function load_and_return_raw_data is defined and loads data from files)
string = "data_config_" * "20" * "ma.jl"
include(string)
FP, time = load_and_return_raw_data(fp_file, ld_file, fp_direct_file, ld_direct_file)

axs[1].set_xlabel("time (s)")
axs[1].set_ylabel("Signal (V)")
axs[1].plot(time[1:100:end], FP[1:100:end], linestyle="-", alpha=1, linewidth=0.4, label="Data Points")

# Add vertical line
inset_time = 0.55
axs[1].axvline(x=inset_time, color="r", linestyle="--",linewidth=0.4)

# Customize tick labels for the first subplot
axs[2].set_xlabel("time (s)")

# Find the index for the inset time and plot a zoomed-in section on the second subplot
idx = findfirst(x -> x > inset_time, time)
axs[2].plot(time[idx:1:idx+200000], FP[idx:1:idx+200000], linestyle="-", alpha=1, linewidth=0.4)
ylim = axs[1].get_ylim()
axs[2].set_ylim(ylim)

# Display the combined plot
plt.tight_layout(pad=0.4)

plt.savefig("./matplotlibfigures/raw_data2",dpi=800)

