include("gas_experiment_functions.jl") 
load_det_data=false

global const sampling_rate=1e6
global const mod_rate=5e3 
println("NB set modulation rate. current = $mod_rate")
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)


# det_data_path="data/20241125/5min_gas.mat"
# ld_data_path="data/20241127/LD_temp_file_with_modulation_10min_gas.txF

###########################3#################################################20241127 bad gas setting_
# det_data_path="data/20241127/mod_10/10min_gas_mod_10.mat" #??
# ld_data_path="data/20241127/mod_10/10min_gas_mod_10.txt"

# det_data_path="data/20241127/2Mhz_mod_20/10min_gas_mod_20_2MHz.mat" #bin exists ?? dividing by time error
# ld_data_path="data/20241127/2Mhz_mod_20/10min_gas_mod_20_2MHz.txt"

##############################################################################################################20241128


# det_data_path="data/20241128/gas_mod10.mat" 
# ld_data_path="data/20241128/gas_mod10.txt"

# det_data_path="data/20241128/gas_mod20.mat"
# ld_data_path="data/20241128/gas_mod20.txt"

#################################################################1129 pressure

# det_data_path="data/20241129/01bar.mat"
# ld_data_path="data/20241129/01bar.txt"


#  det_data_path="data/20241129/02bar.mat" 
#  ld_data_path="data/20241129/02bar.txt"


# det_data_path="data/20241129/03bar.mat" 
# ld_data_path="data/20241129/03bar.txt"
#################################################################20241202 long 01 bar

# det_data_path="data/20241202/01_bar_long.mat" 
# ld_data_path="data/20241202/01bar_long.txt"

# det_data_path="data/20241202/03bar_scan.mat" #long gas flow time >10 min. 20mv 19 pm 2 deg.
# ld_data_path="data/20241202/03bar_scan.txt"

# det_data_path="data/20241202/03bar_stab.mat" #long gas flow time >10 min. 20mv 19 pm 2 deg.

##############################################################################20241211

# det_data_path="data/20241211/20_long.mat" 
# ld_data_path="data/20241211/20_long.txt"


# det_data_path="data/20241211/40_long.mat" 
# ld_data_path="data/20241211/40_long.txt"

##############################################################################20241220 2 and k kHz, with and without loss.

# det_data_path="data/20241220/2khz_20.mat" 
# ld_data_path="data/20241220/2khz_20.txt"

# det_data_path="data/20241220/2khz_20_loss.mat" 
# ld_data_path="data/20241220/2khz_20_with_loss.txt"


det_data_path="data/20241220/5khz_20.mat" 
ld_data_path="data/20241220/5khz_20.txt"

global const mod_rate=5e3 #NB!! for 5khz
# det_data_path="data/20241220/5khz_20_loss.mat" 
# ld_data_path="data/20241220/5khz_20_with_loss.txt"
###########################################################################################################################################################################start
calculate_Lprime_over_L=false 
# # global const mod_rate=5e3 NB!! double check this and set where located.

ld_data=load_ld_data(ld_data_path)

ref_chan_dict, sig_chan_dict,  time_chan_dict,  temp_setpoint_dict,  temp_sensor_dict,  time_ld_dict = process_det_data(det_data_path,ld_data,load_det_data)
 
#initialize dictionaries for data storage and saving
L_dict = Dict{Int, Array{Float64}}()
L_prime_dict = Dict{Int, Array{Float64}}()
trig_dict = Dict{Int, Tuple{Float64, Float64}}()
L_temp_axis = Dict{Int, Array{Float64}}()
Lprime_over_L= Dict{Int, Array{Float64}}()
trig_indices_dict= Dict{Int, Array{Float64}}()

# Loop over all keys in the dictionaries
function get_normalization_coefficient(ref_chan,sig_chan)
    downsample = 100
    # Find reference trigger level
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)

    # Find trigger indices
    trig_indices = find_pulse_trig_points(ref_chan, ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

    # Trim channels
    ref_chan, sig_chan, trig_indices = trim_channels(ref_chan, sig_chan, trig_indices, 2)
    return mean(ref_chan[trig_indices[1]:trig_indices[3]]),mean(sig_chan[trig_indices[1]:trig_indices[3]])
end

# ref_norm,sig_norm=get_normalization_coefficient(ref_chan_dict[1],sig_chan_dict[1]) #outdated?

if calculate_Lprime_over_L==true
    for i in keys(ref_chan_dict)
        println("key=$i")
    
        downsample = 100 #speed up finding the trigger level fit.

        # Find reference trigger level
        ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan_dict[i], downsample)

        # Find trigger indices
        trig_indices = find_pulse_trig_points(ref_chan_dict[i] , ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

        # Trim channels
        ref_chan_trim, sig_chan_trim, time_chan_trim, trig_indices_trim = trim_channels(ref_chan_dict[i] , sig_chan_dict[i],time_chan_dict[i], trig_indices, 2)

        ref_chan_trim=ref_chan_trim
        sig_chan_trim=sig_chan_trim

        # Recalculate the reference trigger level for trimmed and normalized data
        ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan_trim, downsample)
        trig_dict[i]= (ref_trig_level_slope, ref_trig_level_intercept)
        ref_chan_dict[i]= ref_chan_trim 
        sig_chan_dict[i] = sig_chan_trim
        # time_chan_dict[i] = time_chan_trim
        # Run garbage collection
        GC.gc()

        # Calculate L and L_prime for the current scan
        Lprime_over_L[i] = calculate_L_prime_over_L_function(ref_chan_trim, sig_chan_trim, trig_indices_trim, ref_trig_level_slope, ref_trig_level_intercept)
        L_time_axis=[time_chan_trim[i] for i in trig_indices_trim[2:2:end-1]] #the mid point of each scan.
        time_indices=[findmin(abs.(time_ld_dict[i].-x))[2] for x in L_time_axis]
        # time_ld_dict[i][time_chan_indices]
        L_temp_axis[i]=temp_sensor_dict[i][time_indices] 
    end
    for key in keys(L_dict)
        if iseven(key)
            Lprime_over_L[key]=reverse(Lprime_over_L[key])
            L_temp_axis[key]=reverse(L_temp_axis[key])

        end
    end

    save(det_data_path[1:end-3]*"Lprime_over_L.jld2",Dict(string(key) => value for (key, value) in Lprime_over_L))
    save(det_data_path[1:end-3]*"_xaxis.jld2",Dict(string(key) => value for (key, value) in L_temp_axis))
else
    for i in keys(ref_chan_dict)
        println("key=$i")
        # Extract the channels for the current key
        
        downsample = 100
        # Find reference trigger level
        ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan_dict[i], downsample)
        #NB trigger will be off for any hooked shaped scan data, i.e when the ld temp goes through a minimum or maximum. Effect will be greatest at the start of the scan, but will have a small effect on all trigger points, due to wrong trig slope.

        # Find trigger indices
       trig_indices = find_pulse_trig_points(ref_chan_dict[i] , ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)
      
        # Trim channels
        ref_chan_trim, sig_chan_trim, time_chan_trim, trig_indices_trim = trim_channels(ref_chan_dict[i] , sig_chan_dict[i],time_chan_dict[i], trig_indices, 2)
        trig_indices_dict[i]=trig_indices_trim
        #  @infiltrate 
        # Normalize channels
        # ref_chan, sig_chan = normalize_channels(ref_chan, sig_chan, trig_indices)
        ref_chan_trim=ref_chan_trim#/ref_norm
        sig_chan_trim=sig_chan_trim#/sig_norm

        # Recalculate the reference trigger level for trimmed and normalized data
        ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan_trim, downsample)
        trig_dict[i]= (ref_trig_level_slope, ref_trig_level_intercept)
        ref_chan_dict[i]= ref_chan_trim 
        sig_chan_dict[i] = sig_chan_trim
        # time_chan_dict[i] = time_chan_trim
        # Run garbage collection
        GC.gc()

        # Calculate L and L_prime for the current scan
        L_dict[i],L_prime_dict[i],_,_ = calculate_L_L_prime_for_scan(ref_chan_trim, sig_chan_trim, trig_indices_trim, ref_trig_level_slope, ref_trig_level_intercept,sig_norm/ref_norm
        )
        L_time_axis=[time_chan_trim[i] for i in trig_indices_trim[2:2:end-1]] #the mid point of each scan.
        time_indices=[findmin(abs.(time_ld_dict[i].-x))[2] for x in L_time_axis]
        # time_ld_dict[i][time_chan_indices]
        L_temp_axis[i]=temp_sensor_dict[i][time_indices]#temp sensor or temp setpoint?      
    end
    for key in keys(L_dict)
        if iseven(key)
            L_dict[key]=reverse(L_dict[key])  
            L_prime_dict[key]=reverse(L_prime_dict[key])  
            L_temp_axis[key]=reverse(L_temp_axis[key])
        end
    end
    
    # p=plot()
    # for key in keys(L_prime_dict)
    #     plot!(p,L_temp_axis[key],L_prime_dict[key],xlabel="temperature",ylabel="transmittance",label="scan $key")
    #     # plot!(p,L_dict[key])
    # end
    # display(p)
    # savefig(det_data_path[1:end-3]*"_L_prime.png")

    # p=plot()
    # for key in keys(L_prime_dict)
    #     plot!(p,L_temp_axis[key],L_prime_dict[key],xlabel="temperature",ylabel="transmittance",label="scan $key")
    #     # plot!(p,L_dict[key])
    # end
    # display(p)
    # savefig(det_data_path[1:end-3]*"_L_prime.png")

    # averages=101 #NB odd number for center
    # p=plot()
    # for key in keys(L_dict)
    #     xdata=L_temp_axis[key][Int(averages/2+0.5):Int(end-averages/2+0.5)]
    #     ydata=rollmean(L_dict[key],averages)
    #     println(length(xdata))
    #     println(length(ydata))
    #     plot!(p,xdata,ydata,xlabel="temperature",ylabel="transmittance",label="scan $key")
    #     # plot!(p,L_dict[key])
    # end
    # display(p)
    # savefig(det_data_path[1:end-3]*"_L_avg_$averages.png")


    # p=plot()
    # for key in keys(L_dict)
    #     xdata=L_temp_axis[key][Int(averages/2+0.5):Int(end-averages/2+0.5)]
    #     ydata=rollmean(L_prime_dict[key],averages)
    #     println(length(xdata))
    #     println(length(ydata))
    #     plot!(p,xdata,ydata,xlabel="temperature",ylabel="transmittance",label="scan $key")
    #     # plot!(p,L_dict[key])
    # end
    # display(p)
    # savefig(det_data_path[1:end-3]*"_L_prime_avg_$averages.png")

    save(det_data_path[1:end-3]*"_L.jld2",Dict(string(key) => value for (key, value) in L_dict))
    save(det_data_path[1:end-3]*"_L_prime.jld2",Dict(string(key) => value for (key, value) in L_prime_dict))
    save(det_data_path[1:end-3]*"_xaxis.jld2",Dict(string(key) => value for (key, value) in L_temp_axis))
end
# using LsqFit

# itp,_,start_wavelength,stop_wavelength=interpolate_hitran()

# function fit_to_hitran(x,L,itp)
    
#      measurement(c, chi, f, x) = c * (1 .- chi .* itp.(x))

#     # Define your model function f(x)
#     f(x) = exp(-x)  # Example function

#     # Define the fitting model with parameters [c, chi]
#     model(p, x) = measurement(p[1], p[2], f, x)


#     # Initial guess for parameters [c, chi]
#     p0 = [1.0, 0.1]

#     # Perform the fit
#     fit = curve_fit(model, p0, x2, y2)

#     # Extract fitted parameters
#     fitted_params = fit.param
#     return fitted_params
# end

