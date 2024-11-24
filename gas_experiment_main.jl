using MAT,DelimitedFiles,Dates,Plots,Statistics,EasyFit

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

#problem: time varying amplitude and 
function normalize_det_data(sig_chan,ref_chan)
    # mean(ref_chan)
end

function divide_scans(sig_chan,ref_chan)
    return sig_chan[1:66000000],ref_chan[1:66000000]
end

function find_ref_trigger_level(ref_chan,downsample)
    ref_chan=ref_chan[1:downsample:end]
    x=[1:length(ref_chan);]*1.0
    trig_fit=EasyFit.fitlinear(x,ref_chan.*1.0)
    ref_trig_level_slope=trig_fit.a/downsample
    ref_trig_level_intercept=trig_fit.b
    return ref_trig_level_slope,ref_trig_level_intercept
end

function find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope)
    half_periods_of_data = round(Int, length(ref_chan) / data_points_per_half_period)
    trig_indices = Vector{Int}(undef, half_periods_of_data) 
    for i in 1:half_periods_of_data
        local_data=ref_chan[(i - 1) * data_points_per_half_period + 1:(i*data_points_per_half_period)]
        local_trig_level=ref_trig_level_slope*i*data_points_per_half_period+ref_trig_level_intercept #not exact but error should be smaller than data spacing
        trig_index=findmin(abs.(local_data.-local_trig_level))[2]
        trig_tolerance=100
        if abs(local_data[trig_index]-local_trig_level)>trig_tolerance
            error("the difference in the trig level and trig point is greater than $trig_tolerance mv")
        end
        trig_index=trig_index + (i - 1) * data_points_per_half_period #ignore the plus 1 due to the addition.
        trig_indices[i]=trig_index
     end 
    return trig_indices
 
end

function L0(i_p,i_m,Δim)
    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m,Δν_h,Δim)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function average_values(local_data,rise_time,fall_time,trig_indices,index)
    i_plus=       
    i_minus=
    for i in 1:sample_periods
        FP_array[:,i]=FP[indices[i]:Int(indices[i]+T-1)]
    end
    for i in 1:averages
        ts=ts_range[i]
        i_plus[:,i]=FP_array[ts,:].+FP_array[Int(ts+T/2),:]
        i_minus[:,i]=FP_array[ts,:].-FP_array[Int(ts+T/2),:]
        # sample_plus=range(ts,step=T,length=sample_periods) #defines evenly spaced samples. Need to correct for drift in the trigger.
        # sample_minus=sample_plus.+Int(T/2)
        # i_plus[:,i]=FP[sample_plus].+FP[sample_minus]
        # i_minus[:,i]=FP[sample_plus].-FP[sample_minus]
    end
    i_plus=mean(i_plus,dims=2)
    i_minus=mean(i_minus,dims=2)
    return i_plus/I_0,i_minus/I_0
end

function calculate_L_L_prime_for_scan(sig_chan,ref_chan,trig_indices)
    trigger_level(x) = ref_trig_level_slope *x .+ ref_trig_level_intercept
    rise_time=7
    fall_time=2 #NB should be equivalent to be safe #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
    
    L= Vector{Float64}(undef, length(trig_indices[1:2:end]))  
    L_prime=Vector{Float64}(undef, length(trig_indices[1:2:end])) 
    ts_range=
    
    for (index,value) in enumerate(trig_indices[1:2:end]) #the start of every period
        local_data=sig_chan[trig_indices[i]:trig_indices[i+2]]
        I_0=trigger_level(period_number*data_points_per_period)
        high_level=mean(ref_chan[trig_indices[i]+rise_time:trig_indices[i+1]-fall_time])
        low_level=mean(ref_chan[trig_indices[i+1]+fall_time:trig_indices[i+2]-rise_time])
        Δim=high_level-low_level/(2) #I_0 =1

        #FP_modulation_sensitiviy=FP_shift/FP_modulation_amplitude#3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
        # d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second  #NB could also be determined in terms of the min and max wavelength
        Δν_h=-3e10 #1/Δim*(-3e9)#(average_ν-ν_0)#<v>-vo/deltap/p  
        
        i_plus, i_minus = average_values(local_data,rise_time,fall_time,trig_indices,index) #
        L[index] = L0.(i_plus, i_minus, Δim) # Assuming Δim is available
        L_prime[index] = L0_prime.(i_plus, i_minus, Δν_h, Δim) # Assuming Δν_h and Δim are available
    end


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
sig_chan,ref_chan=divide_scans(sig_chan,ref_chan)

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)


trig_indices=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope)
ref_chan,sig_chan,trig_indices=trim_channels(ref_chan,sig_chan,trig_indices)
ref_chan,sig_chan=normalize_channels(ref_chan,sig_chan,trig_indices)
downsample=100
ref_trig_level_slope,ref_trig_level_intercept=find_ref_trigger_level(ref_chan,downsample)

calculate_L_L_prime_for_scan(ref_chan,sig_chan,trig_indices)


# ref_chan=ref_chan[1:10000]
# x_values = 1:length(ref_chan)
# trigger_level = ref_trig_level_slope .* x_values .+ ref_trig_level_intercept
# # Plotting
# plot(x_values, ref_chan, label="Signal (ref_chan)", xlabel="Index", ylabel="Amplitude", title="Signal with Trigger Level and Trigger Points")
# plot!(x_values, trigger_level, label="Trigger Level", color=:green, linewidth=2)
# scatter!(trig_indices[1:40], ref_chan[trig_indices[1:40]], color=:red, marker=:circle, label="Trigger Points")

# plot!(legend=:topright)
# ans=find_pulse_trig_points(ref_chan,ref_trig_level_intercept,ref_trig_level_slope)
# plot(ref_chan[trig_indices])
# plot!(trigger_level[trig_indices])
# plot(ref_chan[trig_indices].-trigger_level[trig_indices])