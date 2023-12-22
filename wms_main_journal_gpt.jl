using DelimitedFiles, CSV, DataFrames, Plots, Statistics, MAT, LsqFit
plotlyjs()
include("wms_main_functions.jl")

# Import the data configuration module for the desired data set
include("data_config_1.jl")
#using .DataConfig1

# Function to load and process data
function load_and_process_data()
    include("data_config_1.jl")

    FP_data = matread(fp_file)
    FP = FP_data["data"]
    time = FP_data["time"]
    LD = matread(ld_file)["data"]
    FP_data = nothing  # Free memory

    # Use experiment parameters from exp_params dictionary
  
    data_points_per_period=Int(0.001/sample_interval)

    ## If data triggered on FP, allign the LD pulses and adjust data accordingly. 
    LD=LD*resistor_value
    LD_mean=mean(LD) #find the mean value to set the 'trigger' point. 
    offset=findfirst(x->x>LD_mean,LD)[1] #middle of the first rise i.e trigger on positive pulse
    LD=LD[offset:end] #adjust the voltage according to the resistor
    FP=FP[offset:end]
    time=time[offset:end]
    println(length(LD))
    #derived values
    rise_time=250#?important for I0 #risetime of intensity modulation. 250?
    I_0=LD_mean# mean(LD[1:data_points_per_period])
    Δim= (mean(LD[rise_time:31250-rise_time])-mean(LD[rise_time+31250:data_points_per_period-rise_time]))/(2*I_0) #account for the resistor difference 10kohm ch1 5kohm chan2

    FP_modulation_sensitiviy=FP_shift/FP_modulation_amplitude#3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
    d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second
    Δν_h=1/Δim*(average_ν-ν_0)#<v>-vo/deltap/p  #NB delta im should be the one used in average measurement, but also this shouldn't change.

    #getting the time
    T=data_points_per_period#period in sample numbers 
    ts=rise_time #defines the number of samples before the value stabilizes
    ts_range=Int.(range(rise_time,step=1,stop=31250-rise_time)) #31250 is the number of data points per half modulation cycle# important in the averaging functions. 
    sample_periods=999 #number of modulation cycles to average over
    averages=Int(length(ts_range)) #number of data points to average over in each modulation cycle
    λ=FP_min_wavelength .+ d_lambda_d_t .* range(0,step=0.001,length=500) #NB 499?
    ν=transpose(c/λ)
    indices=trigger_values(LD,T,sample_periods)

    #############################################################################
    #create the timing and avoid drift in the LD (cant count on an exact number of samples per period)

    ##get the 'direct' or reference data (unmodulated) file names ending in nm for no modulation.
    L_direct=matread(fp_direct_file)["data"]
    I0_direct=matread(ld_direct_file)["data"]
    I0_direct=mean(I0_direct)*resistor_value
    L_direct=L_direct/I0_direct #normalize the direct measurement. NB should use its own reference channel if possible. since current setting may not be exactly I_0
    #LD_direct_mean=matread(f[2])
    L_direct=L_direct[offset:end]
    L_direct_avg=LD_average_values(L_direct,indices,sample_periods,T) 
    L_direct_avg=L_direct_avg[1:500] 

    # Ensure all required variables are passed to functions
    i_plus, i_minus = average_values(FP,indices,sample_periods,averages) # No change needed
    L_direct_avg = LD_average_values(L_direct,indices,sample_periods,T) # No change needed

    FP_array, direct_array, LD_array = chirp_estimate(FP, LD, L_direct) # Pass FP, LD, and L_direct

    # lorentzian_fit() and convert_to_dv() are using global variables, ensure they are available
    L_direct_fit, L_prime_direct = lorentzian_fit(λ,L_direct_average)
    dI_dν = convert_to_dv()

    # Pass necessary variables to L0 and L0_prime functions
    L = L0.(i_plus, i_minus, Δim) # Assuming Δim is available
    L_prime = L0_prime.(i_plus, i_minus, Δν_h, Δim) # Assuming Δν_h and Δim are available

    fwhm_method, fwhm_direct = get_fwhm(λ, L, L_direct_avg)

    chirp = get_chirp_estimate(FP_array, direct_array, LD_array)

    println("FWHM method: $fwhm_method, FWHM direct: $fwhm_direct")
    return λ, ν, L, L_prime, L_direct, dI_dν
end

# Main function
λ,ν,L,L_prime,L_direct, dI_dν= load_and_process_data()

