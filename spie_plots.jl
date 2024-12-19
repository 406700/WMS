## compare L and L_direct.
using JLD2,Plots,RollingFunctions, PyPlot
include("gas_experiment_functions.jl")

total_width_in_inches = 3.5
num_subplots = 1
subplot_width_in_inches = total_width_in_inches / num_subplots
subplot_height_in_inches = total_width_in_inches*3/5
fig, axs = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_height_in_inches))


direct_path="data/20241211/direct_long.mat"
mod_path = "data/20241211/20_long.mat" 
mod_path = "data/20241211/20_long.txt"

# # Load data
mod_L = load(mod_path[1:end-3] * "_L.jld2") 
mod_x = load(mod_path[1:end-3] * "_xaxis.jld2") 

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

#     x = rollx[1:5001:end]
#     y = y ./ I0
#     y = y[1::end]
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
half_w = (window_size - 1) ÷ 2

# Define function to check if a file exists
function check_and_load_rolling_mean(filename)
    if isfile(filename)
        println("Loading rolling mean from disk: $filename")
        return load(filename)["data"]
    else
        return nothing
    end
end

# Define function to compute and save the rolling mean
function compute_and_save_rolling_mean(data, window_size, filename)
    # Compute rolling mean
    result = rollmean(data, window_size)
    
    # Save the result to disk
    @show "Saving rolling mean to disk: $filename"
    save(filename, "data", result)
    
    return result
end

# Parameters
window_size = 5001
half_w = (window_size - 1) ÷ 2

key = "1"
y_filename = "y_rolling_mean_$key.jld2"
I0_filename = "I0_rolling_mean_$key.jld2"

# Check if rolling mean for y is already computed and saved
y = check_and_load_rolling_mean(y_filename)
if y == nothing
    y = compute_and_save_rolling_mean(direct_L[key], window_size, y_filename)
end

# Check if rolling mean for I0 is already computed and saved
I0 = check_and_load_rolling_mean(I0_filename)
if I0 == nothing
    I0 = compute_and_save_rolling_mean(direct_I0[key], window_size, I0_filename)
end


# key="1"
# y = rollmean(direct_L[key], window_size)
# I0 = rollmean(direct_I0[key], window_size)

x = collect(1:length(direct_L[key]))[half_w+1:end-half_w]  
y=y[1:5000:end]./I0[1:5000:end] #NB averaged direct line over averaged I0. correct statistics?
I0=nothing
GC.gc()

xstart=mod_x[key][1] #calibrate the x axis between direct and modulated experiments
xend=mod_x[key][end]
x=collect(range(xstart,stop=xend,length=length(direct_L[key]))) 
x= x[half_w+1:end-half_w]    
x=x[1:5000:end]

#calculate the spline fit and its derivative of the normalized direct measurement
spline=Spline1D(x, y, x[1000:1000:end];  k=3) #use knot vectors with larger spacing
yprime=[derivative(spline,x_value) for x_value in x]
# ans=yprime./y #calculate the normalized derivative i.e dα/dν
Plots.plot(x,yprime/maximum(yprime),label="direct")
# yprime=nothing
# y=nothing
GC.gc()

window_size=500
# Lprime_over_L=load(mod_path[1:end-3] * "Lprime_over_L.jld2") 
# L_prime_over_L=rollmean(Lprime_over_L[key],window_size)
# half_w = (window_size - 1) ÷ 2
# rollx = mod_x[key][half_w+1:end-half_w]
L_prime=load(mod_path[1:end-3] * "_L_prime.jld2") 
L_prime=rollmean(L_prime[key],window_size)
half_w = (window_size - 1) ÷ 2
rollx = mod_x[key][half_w+1:end-half_w]
rollx=rollx[1:window_size:end]
L_prime=L_prime[1:window_size:end]
println("NB manual scaling in the derivative comparison")
Plots.plot!(rollx,-L_prime/(maximum(-L_prime))*1.08,label="modulated",ylabel="derivative of T (a.u.)",xlabel="temperature (C)") #NB minus sign.
Plots.savefig("spie_figures/derivative_comparison.png")

# spline2=Spline1D(rollx, -L_prime, rollx[200:400:end];  k=3)
# Plots.plot(spline2.(rollx))


Plots.plot(x,y,label="direct")
Plots.plot!(x,spline.(x),label="spline")
Plots.scatter!(x[1000:1000:end],ones(length(x[1000:1000:end])),label="nodes")

Plots.savefig("spie_figures/spline_vs_direct")


# ax.legend()
# ax.set_xlabel("temperature")  # Replace with appropriate label
# ax.set_ylabel("scaled derivative, uncalibrated")  # Replace with appropriate label
# plt.tight_layout()

# plt.savefig("spie_figures/derivative_comparison.png",dpi=600)

# x=1:1000
# y=x.^2
# sp1d=Spline1D(x, y;  k=3, bc="nearest", s=3.0)
# Plots.plot(sp1d.(x))


# Plots.plot(sp1d.(x[1:5000:end]))
# Plots.plot!(y[1:5000:end])
# 200/0.2


###########################################################################hitran plot


# Main Code
hitran_path = "data/spectraplot/"
files=readdir(hitran_path)
hitran_path=hitran_path*files[1]
det_data_paths=["data/20241211/20_long.mat"]

p1,p2=compare_hitran_derivative_with_adjusted_scale(hitran_path, det_data_paths)

savefig(p1,"/spie_figures/"*det_data_paths[1][6:end-3]*"hitran_L_prime.png") #nb for vector of paths
savefig(p2,"/spie_figures/"*det_data_paths[1][6:end-3]*"hitran_L.png")

display(p1)

function compare_hitran_derivative_with_adjusted_scale(hitran_path, det_data_paths)
    data, header = readdlm(hitran_path, ',', header = true)
    x = data[:, 1]
    y = data[:, 2]
    # matwrite("hitran.mat",Dict("wavenumber"=>x,"transmittance"=>y))
    p=plot(xlabel="wavenumber",ylabel="uncalibrated derivative (a.u.)")
    p2 = plot(xlabel="wavenumber",ylabel="Transmittance")
    χ_shift = 1

    plot!(p2, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN")
    x,y =numerical_derivative(x,y)
    # matwrite("hitran_derivative.mat",Dict("wavenumber"=>x,"derivative_transmittance"=>y))

    plot!(p, x, (1 .-((1 .-y)* χ_shift)), label = "")
    max_deriv=maximum(y)
    xshift = 0.3
    xstretch = 0.8
    ystretch = 1.0
    yshift = 0#
    for (i, det_data_path) in enumerate(det_data_paths)
        exp_line = load(det_data_path[1:end-3] * "_L_prime.jld2")["1"]
        window_size=9
        half_w = (window_size - 1) ÷ 2
       
        exp_line = rollmean(exp_line, window_size) 
        exp_x= load(det_data_path[1:end-3] * "_xaxis.jld2")["1"][half_w+1:end-half_w]  
        exp_x=map_range.(exp_x,minimum(exp_x),maximum(exp_x),minimum(x),maximum(x))
      
        exp_x = stretch(exp_x, xstretch)
        exp_line = stretch(exp_line, ystretch)
        plot!(p, exp_x.+xshift, reverse(exp_line)/maximum(exp_line)*max_deriv)#, label = "Experimental Line $(i): $(det_data_path)")
            
        exp_line2 = load(det_data_path[1:end-3] * "_L.jld2")["1"] 
        exp_line2 =rollmean(exp_line2, window_size) 

        plot!(p2, exp_x.+xshift, reverse(exp_line2))
    end

    return p,p2
end

function map_range(x, in_min, in_max, out_min, out_max)
    return (x - in_min) / (in_max - in_min) * (out_max - out_min) + out_min
end