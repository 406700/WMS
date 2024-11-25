using MAT,DelimitedFiles,Dates,Plots,Statistics,EasyFit,RollingFunctions

#should rise and fall times be equivalent?

function alternate_vectors(a, b)
    n = length(b)  # Length of the shorter vector `b`
    result = Vector{eltype(a)}()  # Create an empty vector of the appropriate type
        for i in 1:n
        push!(result, a[i])  # Add element from `a`
        push!(result, b[i])  # Add element from `b`
    end
    push!(result, a[n+1])  # Add the last remaining element from `a`
    return result
end
function optimized_findfirst(ref_chan,trig_level)
    threshold = trig_level * maximum(ref_chan)  # Precompute the threshold
    for (i, val) in enumerate(ref_chan)  # Iterate through the array with indices
        if val > threshold
            return i
        end
    end
    return nothing  # Return nothing if no element matches
end

function read_det_data(det_data_path,trig_level)
    matfile=matopen("data/"*det_data_path)   #Use with read, write, close, keys, and haskey.
    sig_chan=read(matfile, "AI_Ch1")
    ref_chan=read(matfile, "AI_Ch0")
    close(matfile)

    start_trigger=optimized_findfirst(ref_chan,trig_level)
    stop_trigger=optimized_findfirst(reverse(ref_chan),trig_level)
    stop_trigger=length(ref_chan)-stop_trigger
    sig_chan=sig_chan[start_trigger:stop_trigger]
    ref_chan=ref_chan[start_trigger:stop_trigger]
    return ref_chan,sig_chan
end


function divide_scans(ref_chan, sig_chan,start_index,stop_index)
    return ref_chan[start_index:stop_index],sig_chan[start_index:stop_index]
end

function find_ref_trigger_level(ref_chan,downsample)
    ref_chan=ref_chan[1:downsample:end]
    x=[1:length(ref_chan);]*1.0
    trig_fit=EasyFit.fitlinear(x,ref_chan.*1.0)
    ref_trig_level_slope=trig_fit.a/downsample
    ref_trig_level_intercept=trig_fit.b
    return ref_trig_level_slope,ref_trig_level_intercept
end

function find_pulse_trig_points(ref_chan, ref_trig_level_intercept, ref_trig_level_slope,data_points_per_half_period)
    trig_indices = Int[]
    length_ref_chan = length(ref_chan)
    trig_tolerance = 100  # Adjust as needed, for unnormalized data.

    first_trig_index=findmin(abs.(ref_chan[1:data_points_per_half_period].-ref_trig_level_intercept))[2]
    push!(trig_indices,first_trig_index)
    # 2. Find subsequent trigger points by searching in a window centered ahead
    window_size = 20  # Number of data points to search in each window
    half_window = window_size/2

    while true
        # Expected index is data_points_per_half_period ahead of the last trigger point
        expected_index = trig_indices[end] + data_points_per_half_period
        if expected_index > length_ref_chan
            break  # Exit if expected index exceeds data length
        end

        # Define the search window centered around the expected index
        window_start = Int(expected_index - half_window)
        window_end = Int(min(expected_index + half_window - 1, length_ref_chan)) #keeps in bounds
        local_indices = window_start:window_end
        local_data = ref_chan[local_indices]

        # Compute the local trigger level at the expected index
        local_trig_level = ref_trig_level_slope * expected_index + ref_trig_level_intercept

        # Find the index in the window where the data crosses the trigger level
        differences = abs.(local_data .- local_trig_level)
        min_difference, min_idx = findmin(differences)

        if min_difference > trig_tolerance
            error("Trigger point not found within tolerance at index $expected_index.")
        end

        # Get the actual index in ref_chan
        trig_index = local_indices[min_idx]
        push!(trig_indices, trig_index)
    end
    if maximum(diff(diff(trig_indices))) >5
        error("finding trig indices: variation in period > 5 sample")
    end
    return trig_indices
end



function L0(i_p,i_m,Δim)
    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m,Δν_h,Δim)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function average_values(ref_chan,rise_time,trig_indices,I_0) #already normalized
    sample_length=minimum(diff(trig_indices)) #find the spacing between trigger points, and use the smallest one
    i_plus=zeros(sample_length)       
    i_minus=zeros(sample_length)

    for i in rise_time:sample_length
        i_plus[i]=ref_chan[i]+ref_chan[i+sample_length]
        i_minus[i]=ref_chan[i]-ref_chan[i+sample_length]
    end
    i_plus=mean(i_plus)
    i_minus=mean(i_minus)
    return i_plus/I_0,i_minus/I_0 #NB not exactly I_0 but accounts for the difference in the intensity over time
end

function calculate_L_L_prime_for_scan(ref_chan, sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)
    trigger_level(x) = ref_trig_level_slope *x .+ ref_trig_level_intercept
    rise_time=7
    fall_time=2 #NB should be equivalent to be safe #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
    
    L= []
    L_prime=[]
    i_plus_all=[]
    i_minus_all=[]
    for i in 1:2:(length(trig_indices)-2) #the start of every period
        local_data=sig_chan[trig_indices[i]:trig_indices[i+2]]
        I_0=trigger_level(trig_indices[i+1])
        high_level=mean(ref_chan[trig_indices[i]+rise_time:trig_indices[i+1]-fall_time])
        low_level=mean(ref_chan[trig_indices[i+1]+fall_time:trig_indices[i+2]-rise_time])
        Δim=high_level-low_level/(2*I_0)

        #FP_modulation_sensitiviy=FP_shift/FP_modulation_amplitude#3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
        # d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second  #NB could also be determined in terms of the min and max wavelength
        Δν_h=-3e10 #1/Δim*(-3e9)#(average_ν-ν_0)#<v>-vo/deltap/p  
        
        i_plus, i_minus = average_values(local_data,rise_time,trig_indices[i:i+2],I_0) #
        push!(L,L0(i_plus, i_minus, Δim))
        push!(L_prime, L0_prime(i_plus, i_minus, Δν_h, Δim) )
        push!(i_minus_all,i_minus)
        push!(i_plus_all,i_plus)
    end
    return L,L_prime,i_plus_all,i_minus_all
end

function prompt_for_integer()
        try
            print("Please enter an integer: ")
            input_str = readline()
            x = parse(Int, input_str)
            return x
        catch
            error("Invalid input. Please enter a valid integer.")
        end

end

function trim_channels(ref_chan,sig_chan,trig_indices)
    foo=2000 #amount of data to plot
    x_values = 1:length(ref_chan[1:foo])
    p=plot()
    plot!(p,x_values, ref_chan[1:foo], label="Signal (ref_chan)", xlabel="Index", ylabel="Amplitude", title="Signal with Trigger Level and Trigger Points")
    scatter!(p,trig_indices[1:5], ref_chan[trig_indices[1:5]], color=:red, marker=:circle, label="Trigger Points")
    # gui(p)
    first_trig_point=2#prompt_for_integer()
    trig_indices=trig_indices[first_trig_point:end]
    if iseven(length(trig_indices)) #make sure it is an odd number of trig indices. 3 per period+2 each additional period.
        trig_indices=trig_indices[1:end-1]
    end
    #now trig_indices_must start at 1
    return ref_chan[trig_indices[1]:trig_indices[end]],sig_chan[trig_indices[1]:trig_indices[end]], (trig_indices.-trig_indices[1].+1)

end

function normalize_channels(ref_chan,sig_chan,trig_indices)
    ref_chan=ref_chan./mean(ref_chan[trig_indices[1]:trig_indices[3]])
    sig_chan=sig_chan./mean(sig_chan[trig_indices[1]:trig_indices[3]])
    return ref_chan,sig_chan
end
function normalize_channels_tdlas(ref_chan,sig_chan,length)
    ref_chan=ref_chan/mean(ref_chan[1:length])
    sig_chan=sig_chan/mean(sig_chan[1:length])

    return ref_chan,sig_chan
end
function divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)
    trigger_level(x) = ref_trig_level_slope *x + ref_trig_level_intercept
    trigger_levels=trigger_level.([1:length(ref_chan);])
    ref_chan=ref_chan./trigger_levels
    sign_chan=sig_chan./trigger_levels
    return ref_chan,sig_chan
end
det_data_path="with_modulation.mat"

ld_data_path="LD_temp_file_with_modulation.txt"
ld_data=readdlm("data/"*ld_data_path)
ld_time=DateTime.(ld_data[:,3], "yyyy-mm-ddTHH:MM:SS.sss")
ld_data=Dict("temp_setpoint"=>ld_data[:,1],"temp_sensor"=>ld_data[:,2],"time"=>ld_time)

findmax(ld_data["temp_setpoint"])
max_indices=findall(x->x==21,ld_data["temp_setpoint"])
findmin(ld_data["temp_setpoint"]) 
min_indices=findall(x->x==17.0,ld_data["temp_setpoint"])
scan_points=ld_data["time"][alternate_vectors(min_indices,max_indices)] #get the time stamps for all the max and min points of the temperature scan. 


ref_chan,sig_chan=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
ref_chan, sig_chan=divide_scans(ref_chan, sig_chan,1,6600000)

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)

downsample=1000
ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)
trig_indices=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope,data_points_per_half_period)
ref_chan,sig_chan,trig_indices=trim_channels(ref_chan,sig_chan,trig_indices)
ref_chan,sig_chan=normalize_channels(ref_chan,sig_chan,trig_indices)
ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)#calculate a new reference level for trimmed and normalized data

L,L_prime,i_p,i_m=calculate_L_L_prime_for_scan(ref_chan,sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)


plot(L)
plot(L_prime)
plot(sig_chan[1:1000:end])
plot(i_p)
plot(i_m)