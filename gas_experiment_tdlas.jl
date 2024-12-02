
include("gas_experiment_functions.jl")

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
# scan_time=200
load_det_data=true
# global const data_points_per_half_period=round(Int,data_points_per_period/2)
# readdir("data")

# det_data_path="data/20241125/5min_gas_1125.mat"
# # # ld_data_path="20241125/LD_temp_file_5min_gas.txt"
# det_data_path="data/20241127/direct/5min_gas_no_mod.mat" #working
# ld_data_path="data/20241127/direct/5min_gas_no_mod.txt"
# det_data_path="data/20241127/reference/no_gas_ref2.mat"#nb stopped after 1 cycle, didn't capture full rise.
# ld_data_path="data/20241127/reference/no_gas_reference2.txt"

#######################################################################20241128

det_data_path="data/20241128/gas_no_mod.mat" #bin exists ?? dividing by time error
ld_data_path="data/20241128/gas_no_mod.txt"


if load_det_data == false
    #load laser diode (ld) temp data
    ld_data=load_ld_data(ld_data_path)

    #load and configure to det data
    ref_chan,sig_chan,time_chan=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
    det_turning_points,ld_turning_points=divide_scans_by_time(time_chan,ld_data["temp_setpoint"],ld_data["time"])
    GC.gc()

    #coordinate the scans 
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

#normalize by the 

for key in 1:3#keys(ref_chan_dict)
    ref_chan=ref_chan_dict[key] 
    sig_chan=sig_chan_dict[key] 
    downsample=100 #amount of downsample to create the linear fit of I0

    ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))
    ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
    ref_chan,sig_chan=divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)

    ref_chan_dict[key] = ref_chan
    sig_chan_dict[key] = sig_chan
end
ref_chan=nothing
sig_chan=nothing
GC.gc()

save(det_data_path[1:end-3]*"_L.jld2",Dict(string(key) => value for (key, value) in sig_chan_dict))


# p=plot()
# for key in 1:3#keys(ref_chan_dict)
#     datas=rollmean(ref_chan_dict[key],100)
#     plot!(p,datas[1:100:end])
#     # plot!(p,L_dict[key])
# end
# display(p)

# savefig(det_data_path[1:end-3]*"sig_chan_minus_ref.png")
# datas=rollmean((sig_chan_dict[1].-ref_chan_dict[1]),100)

# ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)
# trig_indices=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope,data_points_per_half_period)
# ref_chan,sig_chan,trig_indices=trim_channels(ref_chan,sig_chan,trig_indices)
# ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))
# ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
# ref_chan,sig_chan=divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)


# L,L_prime,i_p,i_m=calculate_L_L_prime_for_scan(ref_chan,sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)
# plot(rollmean(ref_chan,10000)[1:100:end])
# plot(rollmean(sig_chan,10000)[1:100:end])

