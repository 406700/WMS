
include("gas_experiment_functions.jl") 
load_det_data=true

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)


# det_data_path="data/20241125/5min_gas.mat"
# ld_data_path="data/20241127/LD_temp_file_with_modulation_10min_gas.txt"


det_data_path="data/20241125/10min_mod2.mat" #working
ld_data_path="data/20241125/LD_temp_file_with_modulation_10min_gas.txt"

# det_data_path="data/20241125/5min_gas_1125.mat"
# ld_data_path="data/20241125/LD_temp_file_5min_gas.txt"


# det_data_path="data/20241127/mod_10/10min_gas_mod_10.mat" #??
# ld_data_path="data/20241127/mod_10/10min_gas_mod_10.txt"

# det_data_path="data/20241127/2Mhz_mod_20/10min_gas_mod_20_2MHz.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241127/2Mhz_mod_20/10min_gas_mod_20_2MHz.txt"

##############################################################################################################20241128

# det_data_path="data/20241128/5min_gas_no_mod.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241128/5min_gas_no_mod.txt"


# det_data_path="data/20241128/gas_mod10.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241128/gas_mod10.txt"

# det_data_path="data/20241128/gas_mod20.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241128/gas_mod20.txt"

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
    
    ref_chan_dict = Dict{Int, Vector{Float64}}()
    sig_chan_dict = Dict{Int, Vector{Float64}}()
    time_chan_dict = Dict{Int, Vector{Float64}}()

    number_scans = length(det_turning_points) - 1

    GC.gc()

    for index in 1:number_scans
        key = index  # Using integers as keys
        range = det_turning_points[index]:det_turning_points[index+1]  # Define the range once
        if isodd(index)
            ref_chan_dict[key] = ref_chan[range]
            sig_chan_dict[key] = sig_chan[range]
            time_chan_dict[key] = time_chan[range]

        else
            ref_chan_dict[key] = reverse!(ref_chan[range])
            sig_chan_dict[key] = reverse!(sig_chan[range])
            time_chan_dict[key] = reverse!(time_chan[range])

        end
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
        if isodd(index)
            temp_setpoint_dict[key] = ld_data["temp_setpoint"][range]
            temp_sensor_dict[key] = ld_data["temp_sensor"][range]
        else
            temp_setpoint_dict[key] = reverse!(ld_data["temp_setpoint"][range])
            temp_sensor_dict[key] = reverse!(ld_data["temp_sensor"][range])
        end
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

# Loop over all keys in the dictionaries
for i in keys(ref_chan_dict)
    println(i)
    # Extract the channels for the current key
    ref_chan = ref_chan_dict[i]
    sig_chan = sig_chan_dict[i]
    time_chan = time_chan_dict[i]

    downsample = 100
    # Find reference trigger level
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)

    # Find trigger indices
    trig_indices = find_pulse_trig_points(ref_chan, ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

    # Trim channels
    ref_chan, sig_chan, trig_indices = trim_channels(ref_chan, sig_chan, trig_indices, 2)

    # Normalize channels
    ref_chan, sig_chan = normalize_channels(ref_chan, sig_chan, trig_indices)

    # Recalculate the reference trigger level for trimmed and normalized data
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)
    ref_chan_dict[i]= ref_chan 
    sig_chan_dict[i] = sig_chan
    time_chan_dict[i] = time_chan
    # Run garbage collection
    GC.gc()

    # Calculate L and L_prime for the current scan
    L_dict[i],L_prime_dict[i],_,_ = calculate_L_L_prime_for_scan(ref_chan, sig_chan, trig_indices, ref_trig_level_slope, ref_trig_level_intercept)

    # Store results in the dictionary
    
end
p=plot()
for key in keys(L_prime_dict)
    plot!(p,L_prime_dict[key])
    # plot!(p,L_dict[key])
end
display(p)
savefig(det_data_path[1:end-3]*"_L_prime.png")

p=plot()
for key in keys(L_prime_dict)
    plot!(p,L_dict[key])
    # plot!(p,L_dict[key])
end
display(p)
savefig(det_data_path[1:end-3]*"_L.png")
# plot(sig_chan[1:100:end])
# plot!(sig_chan[Int(1e6):Int(1e6)+500]) 
# plot(i_p)
# plot(sig_chan[1:1000:end])
# plot(i_p)
# plot(i_m)
