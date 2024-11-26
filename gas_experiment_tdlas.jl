
include("gas_experiment_functions.jl")

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
scan_time=200
# global const data_points_per_half_period=round(Int,data_points_per_period/2)
# readdir("data")
data=readdir("data/20241125/")
# data_list=["gas_start2.mat","20241125/5min_gas.mat","20241125/no_gas.mat"]
det_data_path="20241125/"*data[4]

# ld_data_path="LD_temp_file_with_modulation.txt"
# ld_data=readdlm("data/"*ld_data_path)
# ld_time=DateTime.(ld_data[:,3], "yyyy-mm-ddTHH:MM:SS.sss")
# ld_data=Dict("temp_setpoint"=>ld_data[:,1],"temp_sensor"=>ld_data[:,2],"time"=>ld_time)

# findmax(ld_data["temp_setpoint"])
# max_indices=findall(x->x==21,ld_data["temp_setpoint"])
# findmin(ld_data["temp_setpoint"]) 
# min_indices=findall(x->x==17.0,ld_data["temp_setpoint"])
# scan_points=ld_data["time"][alternate_vectors(min_indices,max_indices)] #get the time stamps for all the max and min points of the temperature scan. 


data,a,b=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
ref_chan=view(data["AI_Ch0"],:)
sig_chan=view(data["AI_Ch1"],:)

turning_points=find_scan_turning_points(ref_chan,scan_time)
ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))

ref_chan_dict = Dict{Int, Vector{Float64}}()
sig_chan_dict = Dict{Int, Vector{Float16}}()
number_scans = length(turning_points) - 1

GC.gc()

for index in 1:number_scans
    key = index  # Using integers as keys
    range = turning_points[index]:turning_points[index+1]  # Define the range once
    if isodd(index)
        ref_chan_dict[key] = ref_chan[range]
        sig_chan_dict[key] = @view sig_chan[range]
    else
        ref_chan_dict[key] = reverse!(@view ref_chan[range])
        sig_chan_dict[key] = reverse!(@view sig_chan[range])
    end
end
GC.gc()
# Create directories to store slices
try readdir("ref_chan_data")
catch
    mkdir("ref_chan_data")
    mkdir("sig_chan_data")
end
# Memory-map each slice in the dictionary
for key in keys(ref_chan_dict)
    # Save ref_chan slice
    open("ref_chan_data/$key.bin", "w") do file
        write(file, ref_chan_dict[key])
    end
    # Map the file to an array
    ref_chan_dict[key] = Mmap.mmap(
        "ref_chan_data/$key.bin",
        Vector{eltype(ref_chan_dict[key])},
        length(ref_chan_dict[key])
    )

    # Save sig_chan slice
    open("sig_chan_data/$key.bin", "w") do file
        write(file, sig_chan_dict[key])
    end
    # Map the file to an array
    sig_chan_dict[key] = Mmap.mmap(
        "sig_chan_data/$key.bin",
        Vector{eltype(sig_chan_dict[key])},
        length(sig_chan_dict[key])
    )
end
fsig_chan=nothing
ref_chan=nothing
GC.gc()
for index in 1:number_scans

    key = index  # Convert index to a string for the key
    ref_chan=ref_chan_dict[key] 
    sig_chan=sig_chan_dict[key] 
    downsample=1000
    ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
    ref_chan,sig_chan=divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)

    ref_chan_dict[key] = ref_chan
    sig_chan_dict[key] = sig_chan
end
ref_chan=nothing
sig_chan=nothing
GC.gc()

p=plot()
for i in 1:number_scans
    plot!(p,rollmean(sig_chan_dict[i][1:1000:end],1000))
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

