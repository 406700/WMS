using DelimitedFiles,FiniteDiff,Interpolations,Plots,RollingFunctions,JLD2

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

exp_line=load(det_data_path[1:end-3]*"_L.jld2")
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
key="1"
direct_path = "data/20241128/gas_no_mod.mat"  # Direct measurement file
direct_line = load(direct_path[1:end-3] * "_L.jld2")[key] 
direct_line = rollmean(direct_line, 100)[1:500:end]  # Smooth and downsample

# Modulation data
mod_paths = Dict(
    30 => "data/20241125/10min_mod2.mat",
    20=> "data/20241128/gas_mod20.mat",
    10 => "data/20241128/gas_mod10.mat"
)
# p=plot()
x_direct=[] 

for mod in [20,10] 
    mod_line = load(mod_paths[mod][1:end-3] * "_L.jld2")
    mod_line = mod_line[key] 
    println(size(mod_line))
    # Create x-axis for direct_line to match mod10_line length
     global x_direct=collect(range(start=1,stop=length(mod_line);length=length(direct_line))) #inside the loop to avoid loading, but uses the last loop iteration. not precise due to varying legnth of mod lines.
    plot!(p,1:length(mod_line), mod_line, label="Modulation_"*string(mod))
end

# Plot direct_line
#since the number of samples are greater, the axis needs to be scaled. should fix everything to use the time channel.


plot!(p,x_direct, direct_line, label="Direct Measurement") #NB

# Pressure data
pres_paths = Dict(
    1 => "data/20241129/01bar.mat",
    2 => "data/20241129/02bar.mat"
)

# Loop over pres_paths to read, process, and plot each dataset
key="1"
for (pressure, path) in pres_paths
    pres_line = load(path[1:end-3] * "_L.jld2")
    pres_line = pres_line[key] 
    plot!(p,1:length(pres_line), pres_line, label="Pressure $(pressure) bar")
end
display(p)

