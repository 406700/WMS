
include("gas_experiment_functions.jl") 
global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)
function load_ld_data(ld_data_path)
    ld_data=readdlm("data/"*ld_data_path)
    ld_time=DateTime.(ld_data[:,3], "yyyy-mm-ddTHH:MM:SS.sss")
    ld_time=convert.(Int,Dates.value.(ld_time.-ld_time[1]))#convert to int milliseconds from start
    temp_set=convert.(Float64,ld_data[:,1])
    temp_sens=convert.(Float64,ld_data[:,2])
    ld_data=Dict("temp_setpoint"=>temp_set,"temp_sensor"=>temp_sens,"time"=>ld_time)
    return ld_data
end


det_data_path="data/20241125/10min_mod2.mat"
ld_data_path="20241125/LD_temp_file_with_modulation_10min_gas.txt"
ld_data=load_ld_data(ld_data_path)

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

i=3
ref_chan=ref_chan_dict[i];
sig_chan=sig_chan_dict[i];
time_chan=time_chan_dict[i];
downsample=100
ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)
trig_indices=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope,data_points_per_half_period)
ref_chan,sig_chan,trig_indices=trim_channels(ref_chan,sig_chan,trig_indices,2)
ref_chan,sig_chan=normalize_channels(ref_chan,sig_chan,trig_indices)
ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
GC.gc()

L,L_prime,i_p,i_m=calculate_L_L_prime_for_scan(ref_chan,sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)

plot(ref_chan[1:100:end])
gui(plt)
push!(p,plt)    
# plot(sig_chan[1:1000:end])
# plot(i_p)
# plot(i_m)
