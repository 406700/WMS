include("gas_experiment_functions.jl") 
load_det_data=true

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)


# det_data_path="data/20241125/5min_gas.mat"
# ld_data_path="data/20241127/LD_temp_file_with_modulation_10min_gas.tx


# det_data_path="data/20241125/10min_mod2.mat" #working
# ld_data_path="data/20241125/LD_temp_file_with_modulation_10min_gas.txt"

# det_data_path="data/20241125/5min_gas_1125.mat"
# ld_data_path="data/20241125/LD_temp_file_5min_gas.txt"

###########################3#################################################20241127 bad gas setting_
# det_data_path="data/20241127/mod_10/10min_gas_mod_10.mat" #??
# ld_data_path="data/20241127/mod_10/10min_gas_mod_10.txt"

# det_data_path="data/20241127/2Mhz_mod_20/10min_gas_mod_20_2MHz.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241127/2Mhz_mod_20/10min_gas_mod_20_2MHz.txt"

##############################################################################################################20241128


# det_data_path="data/20241128/gas_mod10.mat" 
# ld_data_path="data/20241128/gas_mod10.txt"

# det_data_path="data/20241128/gas_mod20.mat"
# ld_data_path="data/20241128/gas_mod20.txt"

#################################################################1129 pressure

# det_data_path="data/20241129/01bar.mat"
# ld_data_path="data/20241129/01bar.txt"


#  det_data_path="data/20241129/02bar.mat" 
#  ld_data_path="data/20241129/02bar.txt"


# det_data_path="data/20241129/03bar.mat" 
# ld_data_path="data/20241129/03bar.txt"
#################################################################20241202 long 01 bar

# det_data_path="data/20241202/01bar_long.mat" 
# ld_data_path="data/20241202/01bar_long.txt"

det_data_path="data/20241202/03bar_scan.mat" #long gas flow time >10 min. 20mv 19 pm 2 deg.
ld_data_path="data/20241202/03bar_scan.txt"

# det_data_path="data/20241202/03bar_stab.mat" #long gas flow time >10 min. 20mv 19 pm 2 deg.

ld_data=load_ld_data(ld_data_path)

if load_det_data == false

    
    #load and configure to det data

    ref_chan,sig_chan,time_chan=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
    GC.gc()
    # ref_chan=view(data["AI_Ch0"],:)
    # sig_chan=view(data["AI_Ch1"],:)
    # time_chan=view(data["AI_Ch0_Xms"],:)

    #coordinate the scans 

    det_turning_points,ld_turning_points=divide_scans_by_time(time_chan,ld_data["temp_setpoint"],ld_data["time"])
    #det turning points are the indices of the times on the adc card where the temperature meets its max and min points. 
    #ld turning points are the indices for the ld computer.

    ref_chan_dict = Dict{Int, Vector{Float64}}()
    sig_chan_dict = Dict{Int, Vector{Float64}}()
    time_chan_dict = Dict{Int, Vector{Float64}}()

    number_scans = length(det_turning_points) - 1

    GC.gc()

    for index in 1:number_scans
        key = index  # Using integers as keys
        range = det_turning_points[index]:det_turning_points[index+1]  # Define the range once
        ref_chan_dict[key] = ref_chan[range]
        sig_chan_dict[key] = sig_chan[range]
        time_chan_dict[key] = time_chan[range]

    end
    ref_chan=nothing
    sig_chan=nothing
    GC.gc()

    temp_setpoint_dict = Dict{Int, Vector{Float64}}()
    temp_sensor_dict = Dict{Int, Vector{Float64}}()
    time_ld_dict = Dict{Int, Vector{Float64}}()

    number_scans = length(det_turning_points) - 1
    for index in 1:number_scans
        key = index  # Using integers as keys
        range = ld_turning_points[index]:ld_turning_points[index+1]  # Define the range once
        temp_setpoint_dict[key] = ld_data["temp_setpoint"][range]
        temp_sensor_dict[key] = ld_data["temp_sensor"][range]
        time_ld_dict[key] = ld_data["time"][range]

    end
    ld_data=nothing
    GC.gc()
 # Save all six dictionaries to disk
 
    open(det_data_path[1:end-3]*"bin", "w") do io
        serialize(io, (
            ref_chan_dict,
            sig_chan_dict,
            time_chan_dict,
            temp_setpoint_dict,
            temp_sensor_dict,
            time_ld_dict
        ))
end
else
    # Load all six dictionaries from disk
    (
        ref_chan_dict,
        sig_chan_dict,
        time_chan_dict,
        temp_setpoint_dict,
        temp_sensor_dict,
        time_ld_dict
    ) = open(det_data_path[1:end-3]*"bin", "r") do io
        deserialize(io)
    end
end


L_dict = Dict{Int, Array{Float64}}()
L_prime_dict = Dict{Int, Array{Float64}}()
trig_dict = Dict{Int, Tuple{Float64, Float64}}()
L_temp_axis = Dict{Int, Array{Float64}}()

# Loop over all keys in the dictionaries
function get_normalization_coefficient(ref_chan,sig_chan)
    

    downsample = 100
    # Find reference trigger level
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)

    # Find trigger indices
    trig_indices = find_pulse_trig_points(ref_chan, ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

    # Trim channels
    ref_chan, sig_chan, trig_indices = trim_channels(ref_chan, sig_chan, trig_indices, 2)
    return mean(ref_chan[trig_indices[1]:trig_indices[3]]),mean(sig_chan[trig_indices[1]:trig_indices[3]])
end

ref_norm,sig_norm=get_normalization_coefficient(ref_chan_dict[1],sig_chan_dict[1])

for i in keys(ref_chan_dict)
    println("key=$i")
    # Extract the channels for the current key
    

    downsample = 100
    # Find reference trigger level
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan_dict[i], downsample)
    #NB trigger will be off for any hooked shaped scan data, i.e when the ld temp goes through a minimum or maximum. Effect will be greatest at the start of the scan, but will have a small effect on all trigger points, due to wrong trig slope.

    # Find trigger indices
    trig_indices = find_pulse_trig_points(ref_chan_dict[i] , ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

    # Trim channels
    ref_chan_trim, sig_chan_trim, time_chan_trim, trig_indices_trim = trim_channels(ref_chan_dict[i] , sig_chan_dict[i],time_chan_dict[i], trig_indices, 2)

    # @infiltrate 
    # Normalize channels
    # ref_chan, sig_chan = normalize_channels(ref_chan, sig_chan, trig_indices)
    ref_chan_trim=ref_chan_trim/ref_norm
    sig_chan_trim=sig_chan_trim/sig_norm

    # Recalculate the reference trigger level for trimmed and normalized data
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan_trim, downsample)
    trig_dict[i]= (ref_trig_level_slope, ref_trig_level_intercept)
    ref_chan_dict[i]= ref_chan_trim 
    sig_chan_dict[i] = sig_chan_trim
    # time_chan_dict[i] = time_chan_trim
    # Run garbage collection
    GC.gc()

    # Calculate L and L_prime for the current scan
    L_dict[i],L_prime_dict[i],_,_ = calculate_L_L_prime_for_scan(ref_chan_trim, sig_chan_trim, trig_indices_trim, ref_trig_level_slope, ref_trig_level_intercept)
    L_time_axis=[time_chan_trim[i] for i in trig_indices_trim[2:2:end-1]] #the mid point of each scan.
    time_indices=[findmin(abs.(time_ld_dict[i].-x))[2] for x in L_time_axis]
    # time_ld_dict[i][time_chan_indices]
    L_temp_axis[i]=temp_sensor_dict[i][time_indices]#temp sensor or temp setpoint?      
end
for key in keys(L_dict)
    if iseven(key)
        L_dict[key]=reverse(L_dict[key])  
        L_prime_dict[key]=reverse(L_prime_dict[key])  
        L_temp_axis[key]=reverse(L_temp_axis[key])
    end
end

# p=plot()
# for key in keys(L_prime_dict)
#     plot!(p,L_temp_axis[key],L_prime_dict[key],xlabel="temperature",ylabel="transmittance",label="scan $key")
#     # plot!(p,L_dict[key])
# end
# display(p)
# savefig(det_data_path[1:end-3]*"_L_prime.png")

# p=plot()
# for key in keys(L_prime_dict)
#     plot!(p,L_temp_axis[key],L_prime_dict[key],xlabel="temperature",ylabel="transmittance",label="scan $key")
#     # plot!(p,L_dict[key])
# end
# display(p)
# savefig(det_data_path[1:end-3]*"_L_prime.png")

# averages=101 #NB odd number for center
# p=plot()
# for key in keys(L_dict)
#     xdata=L_temp_axis[key][Int(averages/2+0.5):Int(end-averages/2+0.5)]
#     ydata=rollmean(L_dict[key],averages)
#     println(length(xdata))
#     println(length(ydata))
#     plot!(p,xdata,ydata,xlabel="temperature",ylabel="transmittance",label="scan $key")
#     # plot!(p,L_dict[key])
# end
# display(p)
# savefig(det_data_path[1:end-3]*"_L_avg_$averages.png")


# p=plot()
# for key in keys(L_dict)
#     xdata=L_temp_axis[key][Int(averages/2+0.5):Int(end-averages/2+0.5)]
#     ydata=rollmean(L_prime_dict[key],averages)
#     println(length(xdata))
#     println(length(ydata))
#     plot!(p,xdata,ydata,xlabel="temperature",ylabel="transmittance",label="scan $key")
#     # plot!(p,L_dict[key])
# end
# display(p)
# savefig(det_data_path[1:end-3]*"_L_prime_avg_$averages.png")

save(det_data_path[1:end-3]*"_L.jld2",Dict(string(key) => value for (key, value) in L_dict))
save(det_data_path[1:end-3]*"_L_prime.jld2",Dict(string(key) => value for (key, value) in L_prime_dict))
save(det_data_path[1:end-3]*"_xaxis.jld2",Dict(string(key) => value for (key, value) in L_temp_axis))


