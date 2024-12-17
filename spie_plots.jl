## compare L and L_direct.
using JLD2,Plots,RollingFunctions, PyPlot
include("gas_experiment_functions.jl")

total_width_in_inches = 3.5
num_subplots = 1
subplot_width_in_inches = total_width_in_inches / num_subplots
subplot_height_in_inches = total_width_in_inches*3/5
fig, axs = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_width_in_inches*2*0.8))


direct_path="data/20241211/direct_long.mat"
mod_path = "data/20241211/20_long.mat" 
mod_path = "data/20241211/20_long.txt"

# # Load data
# mod_L = load(mod_path[1:end-3] * "_L.jld2") 
# mod_x = load(mod_path[1:end-3] * "_xaxis.jld2") 

direct_L = load(direct_path[1:end-3] * "sig_direct.jld2")
direct_I0 = load(direct_path[1:end-3] * "ref_direct.jld2")
direct_x = load(direct_path[1:end-3] * "_L_direct_temp.jld2")

# Create the figure and axes
fig, ax = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))

# Plot data
# for key in "1"#keys(direct_L)
#     window_size = 5001
#     y = rollmean(direct_L[key], window_size)
#     I0 = rollmean(direct_I0[key], window_size)
#     half_w = (window_size - 1) ÷ 2
#     rollx = direct_x[key][half_w+1:end-half_w]  

#     x = rollx[1:100:end]
#     y = y ./ I0
#     y = y[1:100:end]
#     ax.plot(x, y, label="direct")


#     window_size=9
#     y=rollmean(mod_L[key],window_size)
#     half_w = (window_size - 1) ÷ 2
#     rollx = mod_x[key][half_w+1:end-half_w]    
#     ax.plot=(rollx,y,label="20mV modulation")
#     GC.gc()
# end

# # Finalize plot
# ax.legend()
# ax.set_xlabel("temperature")  # Replace with appropriate label
# ax.set_ylabel("transmittance")  # Replace with appropriate label
# plt.tight_layout()
# plt.savefig("/spie_figures/20mvmod_vs_direct.png",dpi=600)


###############################################################################comparision of Lprime_over_L and -1/L( dα dν)
using Dierckx
window_size = 5001
key="1"
y = rollmean(direct_L[key], window_size)

# I0 = rollmean(direct_I0[key], window_size)
# half_w = (window_size - 1) ÷ 2
# rollx = direct_x[key][half_w+1:end-half_w]  
# y=direct_L[key]
x=collect(1:5000:length(y))
y=y[1:5000:end]
spline=Spline1D(x, y, x[1000:1000:end];  k=3)
yprime=[derivative(spline,x_value) for x_value in x]
ans=yprime./y
Plots.plot!(ans*1e6)

Lprime_over_L=load(mod_path[1:end-3] * "Lprime_over_L.jld2") #Lprime_over_L.jld2"Lprime_over_L.jld2"Lprime_over_L.jld2"Lprime_over_L.jld2"
L_prime_over_L=rollmean(Lprime_over_L[key],10)[1:10:end]
Plots.plot!(L_prime_over_L)
# plot(direct_x[key][1:3000:end], direct_L[key][1:3000:end])

# x=1:1000
# y=x.^2
# sp1d=Spline1D(x, y;  k=3, bc="nearest", s=3.0)
# Plots.plot(sp1d.(x))


# Plots.plot(sp1d.(x[1:5000:end]))
# Plots.plot!(y[1:5000:end])
# 200/0.2