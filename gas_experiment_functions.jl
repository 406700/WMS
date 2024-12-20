using MAT,DelimitedFiles,Dates,Plots,Statistics,EasyFit,RollingFunctions,Mmap,Serialization,Infiltrator,RollingFunctions,JLD2,Interpolations
#should rise and fall times be equivalent?
function optimized_findfirst(ref_chan,trig_level)
    threshold = trig_level * maximum(ref_chan)  # Precompute the threshold
    for (i, val) in enumerate(ref_chan)  # Iterate through the array with indices
        if val > threshold
            return i
        end
    end
    return nothing  # Return nothing if no element matches
end
function average_vector_in_chunks(x, y, chunk_size)
    num_chunks = Int(floor(length(y)/chunk_size))
    remainder = length(y) - chunk_size*num_chunks

    if iseven(chunk_size)
        error("even nun-chucks")
    end

    println(remainder)
    if isodd(remainder) && remainder > 1
        half_remainder = (remainder - 1) ÷ 2
        y = y[half_remainder+1:end-half_remainder]
        x = x[half_remainder+1:end-half_remainder]

    elseif remainder > 1
        half_remainder = remainder ÷ 2
        y = y[half_remainder+1:end-half_remainder]
        x = x[half_remainder+1:end-half_remainder]

    else
        y = y[1:end-1]
        x = x[1:end-1]
    end

    avg = [mean(y[(i-1)*chunk_size+1 : i*chunk_size]) for i in 1:num_chunks]
    x_mid = [x[(i-1)*chunk_size + (chunk_size + 1) ÷ 2] for i in 1:num_chunks]

    return x_mid, avg
end
function nearest_odd(x::Float64)
    # Round to the nearest integer
    i = round(Int, x)
    
    # If it's already odd, return it
    if isodd(i)
        return i
    end

    # Otherwise, determine which adjacent odd integer is closer
    lower_odd = i - 1
    upper_odd = i + 1

    if abs(x - lower_odd) < abs(x - upper_odd)
        return lower_odd
    else
        return upper_odd
    end
end
# function read_det_data(det_data_path,trig_level)
#     data=matread("data/"*det_data_path) 
#     delete!(data, "Timestamps_us\0\0\0")
#     delete!(data, "AI_Ch1_Xms")
#     GC.gc()
#     ref_chan=view(data["AI_Ch0"],:)
#     start_trigger=optimized_findfirst(ref_chan,trig_level)
#     stop_trigger=optimized_findfirst(reverse(ref_chan),trig_level)
#     stop_trigger=length(ref_chan)-stop_trigger
#     # for key in keys(data)
#     #     data[key] = data[key][start_trigger:stop_trigger]
#     # end
#     return data,stop_trigger,start_trigger
#     # sig_chan=sig_chan[start_trigger:stop_trigger]
#     # ref_chan=ref_chan[start_trigger:stop_trigger]
#     # return ref_chan,sig_chan
# end
function read_det_data(det_data_path,trig_level)
    
    matfile=matopen(det_data_path)
    ref_chan=read(matfile, "AI_Ch0")[:]
    start_trigger=optimized_findfirst(ref_chan,trig_level)
    stop_trigger=optimized_findfirst(reverse(ref_chan),trig_level)
    stop_trigger=length(ref_chan)-stop_trigger
    GC.gc()
    sig_chan=read(matfile,"AI_Ch1")[start_trigger:stop_trigger]
    time_chan=read(matfile,"AI_Ch0_Xms")[start_trigger:stop_trigger]
    close(matfile)
    return ref_chan[start_trigger:stop_trigger],sig_chan,time_chan.-time_chan[1]

end
function divide_scans_by_time(time_chan,temp_data,temp_time)
    #time must be in same format, not the trigger on the det data sets the accuracy of the allignment.
    #also, any drift in clocks will not be accounted for, unless normalizing the ellapsed time for both data sets. 

    #assumes scan starts at a minima and ends at a minima
    min=minimum(temp_data)
    max=maximum(temp_data)
    max_points=findall(x->x==max,temp_data)
    min_points=findall(x->x==min,temp_data)
    turning_points=vcat(min_points,max_points)
    turning_points=sort(turning_points) # orders the points by time.

    time_points=temp_time[turning_points]
    time_chan_indices=[findmin(abs.(time_chan.-x))[2] for x in time_points]
    return time_chan_indices,turning_points #
end
function find_scan_turning_points(ref_chan,scan_period)
    data_points_per_half_scan_period=Int(scan_period*1/2*sampling_rate) 
    turning_points = Int[]

    first_turning_point=findmin(ref_chan[1:Int(data_points_per_half_scan_period*1.2)])[2]
    push!(turning_points,first_turning_point)
    # 2. Find subsequent trigger points by searching in a window centered ahead
    window_size = 0.1*data_points_per_half_scan_period  # Number of data points to search in each window
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
function load_ld_data(ld_data_path)
    ld_data=readdlm(ld_data_path)
    ld_time=DateTime.(ld_data[:,3], "yyyy-mm-ddTHH:MM:SS.sss")
    ld_time=convert.(Int,Dates.value.(ld_time.-ld_time[1]))#convert to int milliseconds from start
    temp_set=convert.(Float64,ld_data[:,1])
    temp_sens=convert.(Float64,ld_data[:,2])
    ld_data=Dict("temp_setpoint"=>temp_set,"temp_sensor"=>temp_sens,"time"=>ld_time)
    return ld_data
end
# function read_det_data(det_data_path,trig_level)
#     matfile=matopen("data/"*det_data_path)   #Use with read, write, close, keys, and haskey.
#     sig_chan=read(matfile, "AI_Ch1")
#     ref_chan=read(matfile, "AI_Ch0")
#     close(matfile)

#     start_trigger=optimized_findfirst(ref_chan,trig_level)
#     stop_trigger=optimized_findfirst(reverse(ref_chan),trig_level)
#     stop_trigger=length(ref_chan)-stop_trigger
#     sig_chan=sig_chan[start_trigger:stop_trigger]
#     ref_chan=ref_chan[start_trigger:stop_trigger]
#     return ref_chan,sig_chan
# end

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
    trig_tolerance = Int(250*sampling_rate/1e6)  # Adjust as needed, for unnormalized data.

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
            #error("Trigger point not found within tolerance at index $expected_index. sig=$(local_data[min_idx]), trig level=$local_trig_level")
        end

        # Get the actual index in ref_chan
        trig_index = local_indices[min_idx]
        push!(trig_indices, trig_index)
    end
    if maximum(diff(diff(trig_indices))) >(5*sampling_rate/1e6)
        #error("finding trig indices: variation in period $(maximum(diff(diff(trig_indices))))> 5 sample")
    end
    return trig_indices
end



function L0(i_p,i_m,Δim)
    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m,Δν_h,Δim)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function calculate_iplus_iminus(ref_chan,rise_time,trig_indices,I_0) #already normalized
    # @infiltrate
    sample_length=minimum(diff(trig_indices)) #find the spacing between trigger points, and use the smallest one to keep in bounds
    
    i_plus=ref_chan[rise_time:sample_length-rise_time].+ref_chan[rise_time+sample_length:Int(2*sample_length-rise_time)]
    i_minus=ref_chan[rise_time:sample_length-rise_time].-ref_chan[rise_time+sample_length:Int(2*sample_length-rise_time)]
    i_plus=mean(i_plus)
    i_minus=mean(i_minus)
   
    i_plus=mean(i_plus)/I_0
    i_minus=mean(i_minus)/I_0
    return i_plus,i_minus #NB not exactly I_0 but accounts for the difference in the intensity over time
end
function calculate_L_prime_over_L_function(ref_chan, sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept)
    trigger_level(x) = ref_trig_level_slope *x .+ ref_trig_level_intercept
    rise_time=7
    fall_time=7 #NB should be equivalent to be safe #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
    
    Lprime_over_L= []
    
    for i in 1:2:(length(trig_indices)-2) #the start of every period
        
        local_data=sig_chan[trig_indices[i]:trig_indices[i+2]]
        I_0=mean(ref_chan[trig_indices[i]:trig_indices[i+2]] )

        if I_0<0
            error("I_0 less than 1")
        end
        high_level=mean(ref_chan[trig_indices[i]+rise_time:trig_indices[i+1]-fall_time])
        low_level=mean(ref_chan[trig_indices[i+1]+fall_time:trig_indices[i+2]-rise_time])
        if high_level<low_level
            error("high level below low level")
        end

        Δim=(high_level-low_level)/(2*I_0)

        i_plus, i_minus = calculate_iplus_iminus(local_data,rise_time,trig_indices[i:i+2],1.0) #unnormalized
        ans=(-Δim*i_plus+i_minus)/(i_plus-Δim*i_minus)
        push!(Lprime_over_L,ans)
        

    end
    return Lprime_over_L
end
function calculate_L_L_prime_for_scan(ref_chan, sig_chan,trig_indices,ref_trig_level_slope,ref_trig_level_intercept,norm)
    trigger_level(x) = ref_trig_level_slope *x .+ ref_trig_level_intercept
    rise_time=7
    fall_time=7 #NB should be equivalent to be safe #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
    
    L= []
    L_prime=[]
    i_plus_all=[]
    i_minus_all=[]
    I0_all=[]
    for i in 1:2:(length(trig_indices)-2) #the start of every period
        # @infiltrate
        #calculate a new ΔIm and I0
        local_data=sig_chan[trig_indices[i]:trig_indices[i+2]]
        # I_0=trigger_level(trig_indices[i+1])
        I_0=mean(ref_chan[trig_indices[i]:trig_indices[i+2]] )
        if I_0<0
          
            error("I_0 less than 1")
        end
        high_level=mean(ref_chan[trig_indices[i]+rise_time:trig_indices[i+1]-fall_time])
        low_level=mean(ref_chan[trig_indices[i+1]+fall_time:trig_indices[i+2]-rise_time])
        if high_level<low_level
           
            error("high level below low level")
        end
        Δim=(high_level-low_level)/(2*I_0)

        #FP_modulation_sensitiviy=FP_shift/FP_modulation_amplitude#3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
        # d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second  #NB could also be determined in terms of the min and max wavelength
        Δν_h=-1#-3e10 #1/Δim*(-3e9)#(average_ν-ν_0)#<v>-vo/deltap/p  
        
        i_plus, i_minus = calculate_iplus_iminus(local_data,rise_time,trig_indices[i:i+2],I_0*norm) #
        push!(L,L0(i_plus, i_minus, Δim))
        push!(L_prime, L0_prime(i_plus, i_minus, Δν_h, Δim) )
        push!(i_minus_all,i_minus)
        push!(i_plus_all,i_plus)
        push!(I0_all,I_0)

    end
    # I_minus=i_minus_all.*I0_all
    # I_plus=(i_plus_all.*I0_all)
    # matwrite("i_pm",Dict("i_minus"=>i_minus_all,"i_plus"=>i_plus_all))
    # matwrite("I_pm",Dict("I_minus"=>I_minus,"I_plus"=>I_plus))

    return L,L_prime,i_plus_all,i_minus_all
end

function prompt_for_integer()
        try
            print("Please enter an integer: ")
            input_str = readline()
            x = parse(Int, input_str)
            return x
        catch
            #error("Invalid input. Please enter a valid integer.")
        end

end

function trim_channels(ref_chan,sig_chan,trig_indices,first_trigger)
    # foo=2000 #amount of data to plot
    # x_values = 1:length(ref_chan[1:foo])
    # p=plot()
    # plot!(p,x_values, ref_chan[1:foo], label="Signal (ref_chan)", xlabel="Index", ylabel="Amplitude", title="Signal with Trigger Level and Trigger Points")
    # scatter!(p,trig_indices[1:5], ref_chan[trig_indices[1:5]], color=:red, marker=:circle, label="Trigger Points")
    # # gui(p)
    first_trig_point=0
    if ref_chan[trig_indices[1]]<ref_chan[ trig_indices[1]+2 ] #check the sign of the slope at the trigger. Works due to low sampling rate, i.e no change in slope #NB
        first_trig_point=1
    else
        first_trig_point=2
    end
    println(first_trig_point)
    # first_trig_point=first_trigger#prompt_for_integer()
    trig_indices=trig_indices[first_trig_point:end]
    if iseven(length(trig_indices)) #make sure it is an odd number of trig indices. 3 per period+2 each additional period.
        trig_indices=trig_indices[1:end-1]
    end
    #now trig_indices_must start at 1
   ref_chan=ref_chan[trig_indices[1]:trig_indices[end]]
   sig_chan=sig_chan[trig_indices[1]:trig_indices[end]]
   trig_indices=(trig_indices.-trig_indices[1].+1)
   return ref_chan,sig_chan,trig_indices
end
function trim_channels(ref_chan,sig_chan,time_chan,trig_indices,first_trigger)
    # foo=2000 #amount of data to plot
    # x_values = 1:length(ref_chan[1:foo])
    # p=plot()
    # plot!(p,x_values, ref_chan[1:foo], label="Signal (ref_chan)", xlabel="Index", ylabel="Amplitude", title="Signal with Trigger Level and Trigger Points")
    # scatter!(p,trig_indices[1:5], ref_chan[trig_indices[1:5]], color=:red, marker=:circle, label="Trigger Points")
    # # gui(p)
    first_trig_point=0
    if ref_chan[trig_indices[1]]<ref_chan[ trig_indices[1]+2 ] #check the sign of the slope at the trigger. Works due to low sampling rate, i.e no change in slope #NB
        first_trig_point=1
    else
        first_trig_point=2
    end
    println(first_trig_point)
    # first_trig_point=first_trigger#prompt_for_integer()
    trig_indices=trig_indices[first_trig_point:end]
    if iseven(length(trig_indices)) #make sure it is an odd number of trig indices. 3 per period+2 each additional period.
        trig_indices=trig_indices[1:end-1]
    end
    #now trig_indices_must start at 1
   ref_chan=ref_chan[trig_indices[1]:trig_indices[end]]
   sig_chan=sig_chan[trig_indices[1]:trig_indices[end]]
   time_chan=time_chan[trig_indices[1]:trig_indices[end]]

   trig_indices=(trig_indices.-trig_indices[1].+1)
   return ref_chan,sig_chan,time_chan,trig_indices
end

# function normalize_coefficient(ref_chan,sig_chan,trig_indices)
#     # ref_mean=sum(ref_chan[trig_indices[1]:trig_indices[3]])/(trig_indices[3]-trig_indices[1])
#     # sig_mean=sum(sig_chan[trig_indices[1]:trig_indices[3]])/(trig_indices[3]-trig_indices[1])

#     coefficient=mean(sig_chan[trig_indices[1]:trig_indices[3]])/mean(sig_chan[trig_indices[1]:trig_indices[3]])
#     return coefficient
# end
function process_det_data(det_data_path::String, ld_data, load_det_data::Bool)
    if !load_det_data
        # Load and configure to det data
        ref_chan, sig_chan, time_chan = read_det_data(det_data_path, 0.995)  # trigger level
        GC.gc()

        # Coordinate the scans
        det_turning_points, ld_turning_points = divide_scans_by_time(time_chan, ld_data["temp_setpoint"], ld_data["time"])

        # Initialize dictionaries
        ref_chan_dict = Dict{Int, Vector{Float64}}()
        sig_chan_dict = Dict{Int, Vector{Float64}}()
        time_chan_dict = Dict{Int, Vector{Float64}}()

        number_scans = length(det_turning_points) - 1
        GC.gc()

        # Populate ref_chan_dict, sig_chan_dict, and time_chan_dict
        for index in 1:number_scans
            key = index  # Using integers as keys
            range = det_turning_points[index]:det_turning_points[index + 1]
            ref_chan_dict[key] = ref_chan[range]
            sig_chan_dict[key] = sig_chan[range]
            time_chan_dict[key] = time_chan[range]
        end
        ref_chan = nothing
        sig_chan = nothing
        GC.gc()

        # Initialize more dictionaries
        temp_setpoint_dict = Dict{Int, Vector{Float64}}()
        temp_sensor_dict = Dict{Int, Vector{Float64}}()
        time_ld_dict = Dict{Int, Vector{Float64}}()

        # Populate temp_setpoint_dict, temp_sensor_dict, and time_ld_dict
        for index in 1:number_scans
            key = index  # Using integers as keys
            range = ld_turning_points[index]:ld_turning_points[index + 1]
            temp_setpoint_dict[key] = ld_data["temp_setpoint"][range]
            temp_sensor_dict[key] = ld_data["temp_sensor"][range]
            time_ld_dict[key] = ld_data["time"][range]
        end
        ld_data = nothing
        GC.gc()

        # Save all six dictionaries to disk
        save_path = det_data_path[1:end - 3] * "bin"
        open(save_path, "w") do io
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
        save_path = det_data_path[1:end - 3] * "bin"
        (
            ref_chan_dict,
            sig_chan_dict,
            time_chan_dict,
            temp_setpoint_dict,
            temp_sensor_dict,
            time_ld_dict
        ) = open(save_path, "r") do io
            deserialize(io)
        end
    end

    return (
        ref_chan_dict,
        sig_chan_dict,
        time_chan_dict,
        temp_setpoint_dict,
        temp_sensor_dict,
        time_ld_dict
    )
end

function normalize_channels(ref_chan,sig_chan,trig_indices)
    # ref_mean=sum(ref_chan[trig_indices[1]:trig_indices[3]])/(trig_indices[3]-trig_indices[1])
    # sig_mean=sum(sig_chan[trig_indices[1]:trig_indices[3]])/(trig_indices[3]-trig_indices[1])

    ref_chan=ref_chan./mean(sig_chan[trig_indices[1]:trig_indices[3]])
    sig_chan=sig_chan./mean(sig_chan[trig_indices[1]:trig_indices[3]])
    return ref_chan,sig_chan
end
function normalize_channels_tdlas(ref_chan,sig_chan,avg_length)
    ref_chan=ref_chan/mean(ref_chan[1:Int(avg_length)])
    sig_chan=sig_chan/mean(sig_chan[1:Int(avg_length)])

    return ref_chan,sig_chan
end
function divide_scan_by_slope(ref_chan,sig_chan,ref_trig_level_slope,ref_trig_level_intercept)
    trigger_level(x) = ref_trig_level_slope *x + ref_trig_level_intercept
    trigger_levels=trigger_level.([1:length(ref_chan);])
    ref_chan=ref_chan./trigger_levels
    sig_chan=sig_chan./trigger_levels
    return ref_chan,sig_chan
end
function find_scan_turning_points(ref_chan,scan_period)
    data_points_per_half_scan_period=Int(scan_period*1/2*sampling_rate) 
    turning_points = Int[]

    first_turning_point=findmin(ref_chan[1:Int(data_points_per_half_scan_period*1.2)])[2]
    push!(turning_points,first_turning_point)
    # 2. Find subsequent trigger points by searching in a window centered ahead
    window_size = 0.1*data_points_per_half_scan_period  # Number of data points to search in each window
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
function interpolate_hitran()
    hitran_dict=matread("hitran.mat")
    hitran_wavelength=reverse(+1 ./ hitran_dict["wavenumber"]*1e7)
    itp=LinearInterpolation(hitran_wavelength,reverse(hitran_dict["transmittance"]))
    hitran_dict=nothing
    GC.gc()
    start=minimum(hitran_wavelength)
    stop=maximum(hitran_wavelength)

    hitran_dict=matread("hitran_derivative.mat")
    itp2=LinearInterpolation(hitran_wavelength,reverse(hitran_dict["derivative_transmittance"]))

    return itp,itp2,start,stop
end
function numerical_derivative(x,y)

    line_derivative = zeros(Float64, length(y))  
    for i in 2:Int(length(x)-1)
         # Central difference formula
         line_derivative[i]  = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
     end
     return x,line_derivative
 end