
include("gas_experiment_functions.jl")

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
# scan_time=200
load_det_data=true
# global const data_points_per_half_period=round(Int,data_points_per_period/2)
# readdir("data")

# det_data_path="data/20241125/5min_gas_1125.mat"
# # ld_data_path="20241125/LD_temp_file_5min_gas.txt"
det_data_path="data/20241127/direct/5min_gas_no_mod.mat"
ld_data_path="data/20241127/direct/5min_gas_no_mod.txt"
# det_data_path="data/20241127/reference/no_gas_ref2.mat"#nb stopped after 1 cycle, didn't capture full rise.
# ld_data_path="data/20241127/reference/no_gas_reference2.txt"
ld_data=load_ld_data(ld_data_path)

# findmax(ld_data["temp_setpoint"])
# max_indices=findall(x->x==21,ld_data["temp_setpoint"])
# findmin(ld_data["temp_setpoint"]) 
# min_indices=findall(x->x==17.0,ld_data["temp_setpoint"])
# scan_points=ld_data["time"][alternate_vectors(min_indices,max_indices)] #get the time stamps for all the max and min points of the temperature scan. 

if load_det_data == false
    
    #load and configure to det data

    ref_chan,sig_chan,time_chan=read_det_data(det_data_path,0.95) #trigger level, determines when laser powers on and off
    GC.gc()
    # ref_chan=view(data["AI_Ch0"],:)
    # sig_chan=view(data["AI_Ch1"],:)
    # time_chan=view(data["AI_Ch0_Xms"],:)

    #coordinate the scans 

    det_turning_points,ld_turning_points=divide_scans_by_time(time_chan,ld_data["temp_setpoint"],ld_data["time"])
    for key in keys(ld_data)
        ld_data[key]=ld_data[key][ld_turning_points[end-3:end]]
    end
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

# data,a,b=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
# ref_chan=view(data["AI_Ch0"],:)
# sig_chan=view(data["AI_Ch1"],:)

# turning_points=find_scan_turning_points(ref_chan,scan_time)
# ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))

# ref_chan_dict = Dict{Int, Vector{Float64}}()
# sig_chan_dict = Dict{Int, Vector{Float16}}()
# number_scans = length(turning_points) - 1

# GC.gc()

# for index in 1:number_scans
#     key = index  # Using integers as keys
#     range = turning_points[index]:turning_points[index+1]  # Define the range once
#     if isodd(index)
#         ref_chan_dict[key] = ref_chan[range]
#         sig_chan_dict[key] = @view sig_chan[range]
#     else
#         ref_chan_dict[key] = reverse!(@view ref_chan[range])
#         sig_chan_dict[key] = reverse!(@view sig_chan[range])
#     end
# end
# GC.gc()
# # Create directories to store slices
# try readdir("ref_chan_data")
# catch
#     mkdir("ref_chan_data")
#     mkdir("sig_chan_data")
# end
# # Memory-map each slice in the dictionary
# for key in keys(ref_chan_dict)
#     # Save ref_chan slice
#     open("ref_chan_data/$key.bin", "w") do file
#         write(file, ref_chan_dict[key])
#     end
#     # Map the file to an array
#     ref_chan_dict[key] = Mmap.mmap(
#         "ref_chan_data/$key.bin",
#         Vector{eltype(ref_chan_dict[key])},
#         length(ref_chan_dict[key])
#     )

#     # Save sig_chan slice
#     open("sig_chan_data/$key.bin", "w") do file
#         write(file, sig_chan_dict[key])
#     end
#     # Map the file to an array
#     sig_chan_dict[key] = Mmap.mmap(
#         "sig_chan_data/$key.bin",
#         Vector{eltype(sig_chan_dict[key])},
#         length(sig_chan_dict[key])
#     )
# end
# fsig_chan=nothing
# ref_chan=nothing
# GC.gc()

for key in 1:3
    # @infiltrate
    ref_chan=ref_chan_dict[key] 
    sig_chan=sig_chan_dict[key] 
    downsample=1000
    ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))
    ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
    ref_chan,sig_chan=divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)

    ref_chan_dict[key] = ref_chan
    sig_chan_dict[key] = sig_chan
end
ref_chan=nothing
sig_chan=nothing
GC.gc()

p=plot()
for i in 1:3
    plot!(p,(ref_chan_dict[i][1:1000:end]))
end
display(p)
# ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)
# trig_indices=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope,data_points_per_half_period)
# ref_chan,sig_chan,trig_indices=trim_channels(ref_chan,sig_chan,trig_indices)
# ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))
# ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
# ref_chan,sig_chan=divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)


# L,L_prime,i_p,i_m=calculate_L_L_prime_for_scan(ref_chan,sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)
# plot(rollmean(ref_chan,10000)[1:100:end])
# plot(rollmean(sig_chan,10000)[1:100:end])

