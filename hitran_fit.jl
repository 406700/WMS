using DelimitedFiles,FiniteDiff,Interpolations,Plots,RollingFunctions

hitran_path="data/spectraplot/"
#load a line
# det_data_path="data/20241125/10min_mod2.mat" #working
# det_data_path="data/20241128/gas_mod10.mat" #bin exists ?? dividing by time error
det_data_path="data/20241128/gas_no_mod.mat" #direct measurement
direct_measurement=true

files=readdir(hitran_path)
data,header=readdlm(hitran_path*files[1],',',header=true)
x=data[:,1]
y=data[:,2]

shift(x,shift)=x+shift

function stretch(x_axis, scaling_factor)
    x_center = (maximum(x_axis) + minimum(x_axis)) / 2
    return x_center .+ scaling_factor .* (x_axis .- x_center)
end
# interpolated_line=LinearInterpolation(x, y, extrapolation_bc=NaN)
# f(x)=interpolated_line(x)
# # Compute the derivative of y with respect to x
# line_derivative=FiniteDiff.finite_difference_derivative(f(x),x)

exp_line=readdlm("det_data_path[1:end-3]"*"_L")
# exp_line=Vector(exp_line[:,1])
# if direct_measurement==true
#     exp_line=rollmean(exp_line,100)[1:500:end] #500=same as 
# end


# GC.gc()
# exp_x=[range(x[1],stop=x[end],length=length(exp_line));]
# xshift=+0.1
# xstretch=0.7
# ystretch=1
# yshift=-0.01
# exp_x=stretch(exp_x,xstretch)
# exp_x=shift.(exp_x,xshift)
# exp_line2=stretch(exp_line,ystretch)
# exp_line2=shift.(exp_line,yshift)
# χ_shift=0.5
# plot(x,y*χ_shift.+.49)
# plot!(exp_x,(exp_line2))


###############################################comparing the direct and derived lines 

# direct_path="data/20241128/gas_no_mod.mat" #direct measurement
# direct_line=readdlm(direct_path[1:end-3]*"_L")
# direct_line=Vector(direct_line[:,1])
# direct_line=rollmean(direct_line,100)[1:500:end] #500=same as 

# mod_paths=Dict(30=>"data/20241125/10min_mod2.mat",10=>"data/20241128/gas_mod10.mat")
# pres_paths=Dict(01=>"data/20241129/01bar.mat",02=>"data/20241129/02bar.mat")
# mod10_line=readdlm(mod_paths[10][1:end-3]*"_L")
# mod10_line=Vector(mod10_line[:,1])


# direct_x=[range(1,stop=length(mod10_line);length=length(direct_line));]

# plot(1:length(mod10_line),mod10_line)
# plot!(direct_x,direct_line,label="direct")

# size(mod10_line)

# using DelimitedFiles
# using StatsBase
# using Plots

# Direct measurement
direct_path = "data/20241128/gas_no_mod.mat"  # Direct measurement file
direct_line = readdlm(direct_path[1:end-3] * "_L")
direct_line = Vector(direct_line[:, 1])
direct_line = rollmean(direct_line, 100)[1:500:end]  # Smooth and downsample

# Modulation data
mod_paths = Dict(
    30 => "data/20241125/10min_mod2.mat",
    10 => "data/20241128/gas_mod10.mat"
)
mod10_line = readdlm(mod_paths[10][1:end-3] * "_L")
mod10_line = Vector(mod10_line[:, 1])

# Create x-axis for direct_line to match mod10_line length
direct_x = range(1, stop=length(mod10_line), length=length(direct_line))

# Initialize the plot with mod10_line
p=plot()
plot!(p,1:length(mod10_line), mod10_line, label="Modulation 10")

# Plot direct_line
plot!(p,direct_x, direct_line, label="Direct Measurement")

# Pressure data
pres_paths = Dict(
    1 => "data/20241129/01bar.mat",
    2 => "data/20241129/02bar.mat"
)

# Loop over pres_paths to read, process, and plot each dataset
for (pressure, path) in pres_paths
    pres_line = readdlm(path[1:end-3] * "_L")
    pres_line = Vector(pres_line[:, 1])
    # pres_line = rollmean(pres_line, 100)[1:500:end]  # Smooth and downsample
    # Create x-axis for pres_line to match mod10_line length
    pres_x = range(1, stop=length(mod10_line), length=length(pres_line))
    # Plot pres_line
    plot!(p,pres_x, pres_line, label="Pressure $(pressure) bar")
end

# Display the plot
display(p)

# Print the size of mod10_line
println("Size of mod10_line: ", size(mod10_line))
