
include("gas_experiment_functions.jl")

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
# scan_time=200
load_det_data=false
# global const data_points_per_half_period=round(Int,data_points_per_period/2)
# readdir("data")

# det_data_path="data/20241125/5min_gas_1125.mat"
# # # ld_data_path="20241125/LD_temp_file_5min_gas.txt"
# det_data_path="data/20241127/direct/5min_gas_no_mod.mat" #working
# ld_data_path="data/20241127/direct/5min_gas_no_mod.txt"
# det_data_path="data/20241127/reference/no_gas_ref2.mat"#nb stopped after 1 cycle, didn't capture full rise.
# ld_data_path="data/20241127/reference/no_gas_reference2.txt"

#######################################################################20241128

# det_data_path="data/20241128/gas_no_mod.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241128/gas_no_mod.txt"

#########################################################################20241211
det_data_path="data/20241211/direct_long.mat"
ld_data_path="data/20241211/direct_long.txt"

# det_data_path="data/20241211/dierct_13.mat" 
# ld_data_path="data/20241211/direct_13_deg.txt"

# det_data_path="data/20241211/no_gas_long.mat" 
# ld_data_path="data/20241211/no_gas_long.txt"

ld_data=load_ld_data(ld_data_path)

ref_chan_dict, sig_chan_dict,  time_chan_dict,  temp_setpoint_dict,  temp_sensor_dict,  time_ld_dict = process_det_data(det_data_path,ld_data,load_det_data)
GC.gc()
function get_normalization_coefficient(ref_chan,sig_chan,data_points_per_period)
   ref_norm=mean(ref_chan[1:data_points_per_period])
   sig_norm=mean(sig_chan[1:data_points_per_period])
   return ref_chan/ref_norm, sig_chan/sig_norm     
end

# if load_det_data == false
#     #load laser diode (ld) temp data
#     ld_data=load_ld_data(ld_data_path)

#     #load and configure to det data
#     ref_chan,sig_chan,time_chan=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
#     det_turning_points,ld_turning_points=divide_scans_by_time(time_chan,ld_data["temp_setpoint"],ld_data["time"])
#     GC.gc()

#     #coordinate the scans 
#     ref_chan_dict = Dict{Int, Vector{Float64}}()
#     sig_chan_dict = Dict{Int, Vector{Float64}}()
#     time_chan_dict = Dict{Int, Vector{Float64}}()

#     number_scans = length(det_turning_points) - 1

#     GC.gc()

#     for index in 1:number_scans
#         key = index  # Using integers as keys
#         range = det_turning_points[index]:det_turning_points[index+1]  # Define the range once
#         if isodd(index)
#             ref_chan_dict[key] = ref_chan[range]
#             sig_chan_dict[key] = sig_chan[range]
#             time_chan_dict[key] = time_chan[range]

#         else
#             ref_chan_dict[key] = reverse!(ref_chan[range])
#             sig_chan_dict[key] = reverse!(sig_chan[range])
#             time_chan_dict[key] = reverse!(time_chan[range])

#         end
#     end
#     ref_chan=nothing
#     sig_chan=nothing
#     GC.gc()


#     temp_setpoint_dict = Dict{Int, Vector{Float64}}()
#     temp_sensor_dict = Dict{Int, Vector{Float64}}()
#     time_ld_dict = Dict{Int, Vector{Float64}}()

#     number_scans = length(det_turning_points) - 1
#     for index in 1:number_scans
#         key = index  # Using integers as keys
#         range = ld_turning_points[index]:ld_turning_points[index+1]  # Define the range once
#         if isodd(index)
#             temp_setpoint_dict[key] = ld_data["temp_setpoint"][range]
#             temp_sensor_dict[key] = ld_data["temp_sensor"][range]
#         else
#             temp_setpoint_dict[key] = reverse!(ld_data["temp_setpoint"][range])
#             temp_sensor_dict[key] = reverse!(ld_data["temp_sensor"][range])
#         end
#     end
#     ld_data=nothing
#     GC.gc()
#  # Save all six dictionaries to disk
#     open(det_data_path[1:end-3]*"bin", "w") do io
#         serialize(io, (
#             ref_chan_dict,
#             sig_chan_dict,
#             time_chan_dict,
#             temp_setpoint_dict,
#             temp_sensor_dict,
#             time_ld_dict
#         ))
# end
# else
#     # Load all six dictionaries from disk
#     (
#         ref_chan_dict,
#         sig_chan_dict,
#         time_chan_dict,
#         temp_setpoint_dict,
#         temp_sensor_dict,
#         time_ld_dict
#     ) = open(det_data_path[1:end-3]*"bin", "r") do io
#         deserialize(io)
#     end
# end

L_dict = Dict{Int, Array{Float64}}()
L_temp_axis = Dict{Int, Array{Float64}}()

function nearest_indices_sorted(time_array::AbstractVector{<:Real},   query_times::AbstractVector{<:Real})::Vector{Int}
    # This function assumes both time_array and query_times are sorted in ascending order.
    # Returns a vector of Int indices the same length as query_times.
    # For each element in query_times, it finds the nearest element in time_array.
    n = length(time_array)
    m = length(query_times)
    indices = similar(query_times, Int)

    # If either array is empty, handle gracefully:
    if n == 0
    error("time_array is empty.")
    end
    if m == 0
    return indices # empty vector
    end

    j = 1
    for i in 1:m
    t = query_times[i]

    # Advance j until time_array[j] >= t or we reach the end
    while j < n && time_array[j] < t
    j += 1
    end

    # Now we have a few cases:
    if j == 1
        # t is less than or equal to the first element in time_array
        indices[i] = 1
        elseif j == n && time_array[j] < t
        # t is greater than all elements in time_array
        indices[i] = n
        else
        # Here time_array[j] >= t and j > 1
        # Compare time_array[j] with time_array[j-1] to find the closest
        leftval = time_array[j-1]
        rightval = time_array[j]
        if abs(rightval - t) < abs(leftval - t)
        indices[i] = j
        else
        indices[i] = j - 1
        end
        end
    end
    return indices
end

for key in keys(ref_chan_dict)
    ref_chan=ref_chan_dict[key] 
    sig_chan=sig_chan_dict[key] 
    time_chan=time_chan_dict[key]
    downsample=100 #amount of downsample to create the linear fit of I0

    ref_chan,sig_chan=get_normalization_coefficient(ref_chan,sig_chan,Int(data_points_per_period))
    ref_chan_dict[key] = ref_chan
    sig_chan_dict[key] = sig_chan
    # L_dict[key]=sig_chanref_chan

    L_time_axis=@view(time_chan[1:end])
    time_indices=nearest_indices_sorted(time_ld_dict[key],L_time_axis)
    L_temp_axis[key]=temp_sensor_dict[key][time_indices]#temp sensor or temp setpoint?      
    
end
GC.gc()

for key in keys(L_dict)
    if iseven(key)
        L_dict[key]=reverse(L_dict[key])  
        L_temp_axis[key]=reverse(L_temp_axis[key])
    end
end

GC.gc()
p=plot()
for key in keys(L_dict)
    plot!(p,L_temp_axis[key][1:100:end],L_dict[key][1:100:end])
end
display(p)


save(det_data_path[1:end-3]*"sig_direct.jld2",Dict(string(key) => value for (key, value) in sig_chan_dict))
save(det_data_path[1:end-3]*"ref_direct.jld2",Dict(string(key) => value for (key, value) in ref_chan_dict))
save(det_data_path[1:end-3]*"_L_direct_temp.jld2",Dict(string(key) => value for (key, value) in L_temp_axis))
# foo=load((det_data_path[1:end-3]*"_L.jld2"))
# plot(rollmean(foo["1"],1000)[1:10:end])
# det_data_path="data/20241211/20_long.mat" #long gas flow time >10 min. 20mv 19 pm 2 deg.
# foo=load((det_data_path[1:end-3]*"_L.jld2"))
# plot!(foo["1"])

######################################################################################compare direct and modulation temp location
# # Load the datasets
# det_data_path = "data/20241211/direct_long.mat"
# foo1 = load(det_data_path[1:end-3] * "_L.jld2")

# det_data_path = "data/20241211/20_long.mat"
# foo2 = load(det_data_path[1:end-3] * "_L.jld2")
# averages=5000
# # Process the first dataset (downsample and smooth)
# data1 = rollmean(foo1["2"], averages)[1:100:end]

# # Define the common x-range (assuming 0 to 1 for simplicity)
# x_common = range(0, stop=1, length=length(data1))

# # Scale the second dataset's x-axis to the common range

# data2 = rollmean(foo2["2"],Int(floor(averages/500)))
# x2_scaled = range(0, stop=1, length=length(data2))

# # Plot both datasets
# plot(x_common, data1, label="direct", lw=2)  # First dataset
# plot!(x2_scaled, data2, label="modulated", lw=2) 
# savefig(det_data_path[1:end-3]*"direct_comparison2_$averages.png")
# # p=plot()
# # for key in 1:3#keys(ref_chan_dict)
# #     datas=rollmean(ref_chan_dict[key],100)
# #     plot!(p,datas[1:100:end])
# #     # plot!(p,L_dict[key])
# # end
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

