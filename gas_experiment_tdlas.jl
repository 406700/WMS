function find_scan_turning_points(ref_chan)
    data_points_per_half_scan_period=Int(65e6)
    turning_points = Int[]

    first_turning_point=findmin(ref_chan[1:data_points_per_half_scan_period])[2]
    push!(turning_points,first_turning_point)
    # 2. Find subsequent trigger points by searching in a window centered ahead
    window_size = Int(1e7)  # Number of data points to search in each window
    half_window = window_size÷2
    i=1
    while true
        # Expected index is data_points_per_half_period ahead of the last trigger point
        expected_index = turning_points[end] + data_points_per_half_scan_period
        if expected_index > length(ref_chan)
            break  # Exit if expected index exceeds data length
        end

        # Define the search window centered around the expected index
        window_start = Int(expected_index - half_window)
        window_end = Int(min(expected_index + half_window, length(ref_chan))) #keeps in bounds
        local_indices = window_start:window_end
        local_data = ref_chan[local_indices]
        if isodd(i)
            min_difference, min_idx = findmax(local_data)
        else
            min_difference, min_idx = findmin(local_data)
        end
      

    
        # Get the actual index in ref_chan
        turning_point = local_indices[min_idx]
        push!(turning_points, turning_point)
        i=i+1
    end
   
    return turning_points
end

include("gas_experiment_main.jl")

# readdir("data")
# det_data_path="gas_start2.mat"

# ld_data_path="LD_temp_file_with_modulation.txt"
# ld_data=readdlm("data/"*ld_data_path)
# ld_time=DateTime.(ld_data[:,3], "yyyy-mm-ddTHH:MM:SS.sss")
# ld_data=Dict("temp_setpoint"=>ld_data[:,1],"temp_sensor"=>ld_data[:,2],"time"=>ld_time)

# findmax(ld_data["temp_setpoint"])
# max_indices=findall(x->x==21,ld_data["temp_setpoint"])
# findmin(ld_data["temp_setpoint"]) 
# min_indices=findall(x->x==17.0,ld_data["temp_setpoint"])
# scan_points=ld_data["time"][alternate_vectors(min_indices,max_indices)] #get the time stamps for all the max and min points of the temperature scan. 


ref_chan,sig_chan=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
# ref_chan, sig_chan=divide_scans(ref_chan, sig_chan,131900000,(119000000+66000000))

# ref_chan, sig_chan=divide_scans(ref_chan, sig_chan,9000000,length(ref_chan))
turning_points=find_scan_turning_points(ref_chan)

plot(ref_chan[1:100:end])

downsample=1000
# ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)
# trig_indices=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope,data_points_per_half_period)
# ref_chan,sig_chan,trig_indices=trim_channels(ref_chan,sig_chan,trig_indices)
ref_chan,sig_chan=normalize_channels_tdlas(ref_chan,sig_chan,Int(data_points_per_period))
ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data
ref_chan,sig_chan=divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)


# L,L_prime,i_p,i_m=calculate_L_L_prime_for_scan(ref_chan,sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)
plot(rollmean(ref_chan,10000)[1:100:end])
plot(rollmean(sig_chan,10000)[1:100:end])

