
using DelimitedFiles, CSV, DataFrames, Statistics, MAT, LsqFit,Serialization, RollingFunctions,Infiltrator,EasyFit,Interpolations,Roots,LaTeXStrings,Measures,PyPlot,LaTeXStrings #Plots



import PyPlot; const plt = PyPlot;
using PyCall
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
function set_default_font(font_path::String)
    fm = pyimport("matplotlib.font_manager")
    font_prop = fm.FontProperties(fname=font_path)
    font_name = font_prop.get_name()
    rcParams["font.family"] = font_name
end
set_default_font("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")

# Set default figure size
total_width_in_inches = 3.5
aspect_ratio=3/4
# rcParams["text.usetex"]=true #keep off for now due to font uses. only use if necessary.
rcParams["figure.figsize"] = (total_width_in_inches, total_width_in_inches*aspect_ratio)
rcParams["figure.dpi"] = 800
rcParams["lines.linewidth"]=0.4
rcParams["axes.titlesize"]=8
rcParams["xtick.labelsize"]=6
rcParams["ytick.labelsize"]=6
rcParams["legend.fontsize"]=6
rcParams["figure.titlesize"]=8
# rcParams["lines.marker"] = "."  # Set default marker
rcParams["lines.markersize"] = 4 
# plt.rcParams.update({
#     "font.size': 8,
#     'axes.titlesize': 8,
#     'axes.labelsize': 8,
#     'xtick.labelsize': 8,
#     'ytick.labelsize': 8,
#     'legend.fontsize': 8,
#     'figure.titlesize': 8
# })
# gr()
# upscale = 1 #8x upscaling in resolution
# fntsm = Plots.font("Helvetica", pointsize=round(12.0*upscale))
# fntlg = Plots.font("Helvetica", pointsize=round(16.0*upscale))
# default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
# default(size=(width_in_pixels,)) #Plot canvas size


# Function to load and process data
# function load_and_process_data_old(fp_file,ld_file,fp_direct_file,ld_direct_file)
#     FP_data = matread(fp_file)
#     FP = FP_data["data"]
#     time = FP_data["time"]
#     LD = matread(ld_file)["data"]
#     FP_data = nothing  # Free memory

#     # Use experiment parameters from exp_params dictionary
  
#     data_points_per_period=Int(0.001/sample_interval)
#     T=data_points_per_period
#     ## If data triggered on FP, allign the LD pulses and adjust data accordingly. 
#     LD=LD*resistor_value
#     LD_mean=mean(LD) #find the mean value to set the 'trigger' point. 
#     offset_indices= offset_trigger_values(LD,T,2)
#     offset=find_rising_trigger(LD[1:100000],offset_indices,5)#findfirst(x->x>LD_mean,LD)[1] #middle of the first rise i.e trigger on positive pulse
     
#     LD=LD[offset:end] 
#     FP=FP[offset:end]
#     time=time[offset:end]
#     #derived values
#     rise_time=250 #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
#     I_0=LD_mean
#     half_period=31250
#     Δim= (mean(LD[rise_time:half_period-rise_time])-mean(LD[rise_time+half_period:data_points_per_period-rise_time]))/(2*I_0) #account for the resistor difference 10kohm ch1 5kohm chan2

#     Δν_h=1/Δim*(average_ν-ν_0)#<v>-vo/deltap/p  #NB delta im should be the one used in average measurement, but also this shouldn't change.

#     #getting the time
#     T=data_points_per_period#period in sample numbers 
#     ts=rise_time #defines the number of samples before the value stabilizes
#     ts_range=Int.(range(rise_time,step=1,stop=31250-rise_time)) #31250 is the number of data points per half modulation cycle# important in the averaging functions. 
#     sample_periods=999 #number of modulation cycles 
#     averages=Int(length(ts_range)) #number of data points to average over in each modulation cycle
  
#     d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second
#     FP_min_wavelength=c/ν_0-FP_shift/2
#     λ=FP_min_wavelength .+ d_lambda_d_t .* range(0,step=0.001,length=sample_periods)

#     ν=transpose([c]./[λ;])
  
#     indices=trigger_values(LD,T,sample_periods)

#     #############################################################################
#     #create the timing and avoid drift in the LD (cant count on an exact number of samples per period)

#     ##get the 'direct' or reference data / no modulation
#     L_direct=matread(fp_direct_file)["data"]
#     I0_direct=matread(ld_direct_file)["data"]
#     I0_direct=mean(I0_direct)*resistor_value
#     L_direct=L_direct/I0_direct #normalize the direct measurement. NB should use its own reference channel if possible. since current setting may not be exactly I_0

#     L_direct=L_direct[offset:end]  #there is no reason for the offset rather than to have the data the same length as the modulated data. 
#     L_direct_avg=LD_average_values(L_direct,indices,sample_periods,T) 
#     L_direct_avg=L_direct_avg[1:sample_periods] 

#     # Ensure all required variables are passed to functions
#     i_plus, i_minus = average_values(FP,indices,sample_periods,averages,T,ts_range,I_0) # No change needed
#     FP_array, LD_array = make_arrays_old(FP,LD,sample_periods,T,indices)
    
#     #chirp = get_chirp_estimate(FP_array, LD_array,L_direct_avg,λ, M)
#     #p=plot()
#     #plot(p,chirp*1e-9,title="chirp in nm")
#     #display(p)
#     L_direct_fit, L_prime_direct = lorentzian_fit(λ,L_direct_avg)
#     L_prime_direct = convert_to_dv(L_prime_direct,λ)

#     # Pass necessary variables to L0 and L0_prime functions
#     L = L0.(i_plus, i_minus, Δim) # Assuming Δim is available
#     L_prime = L0_prime.(i_plus, i_minus, Δν_h, Δim) # Assuming Δν_h and Δim are available

#     #confirmation of fwhm
#     fwhm_method, fwhm_direct = get_fwhm(λ, L, L_direct_avg)
#     println("FWHM method: $fwhm_method, FWHM direct: $fwhm_direct")

#     return λ, ν, L, L_prime, L_direct_avg, L_prime_direct,FP_array, LD_array#,chirp
# end;
#calculates arraywith the new make_arrays() functions.
function load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file)
    FP_data = matread(fp_file)
    FP = FP_data["data"]
    time = FP_data["time"]
    LD = matread(ld_file)["data"]
    FP_data = nothing  # Free memory

    # Use experiment parameters from exp_params dictionary
    rise_time=250 #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
    data_points_per_period=Int(0.001/sample_interval)
    T=data_points_per_period
    sample_periods=999
    ## If data triggered on FP, allign the LD pulses and adjust data accordingly. 
    LD=LD*resistor_value
    LD_mean=mean(LD) #find the mean value to set the 'trigger' point. 
    if LD_mean ==Inf
        error("infinite error")
    end
    triggers=Int.(trigger(LD,T,sample_periods))
    offset=triggers[1]
    LD=LD[offset:end] 
    FP=FP[offset:end]
    time=time[offset:end]
    #derived values
    
    I_0=LD_mean
    half_period=31250
    Δim= (mean(LD[2*rise_time:half_period-2*rise_time])-mean(LD[2*rise_time+half_period:data_points_per_period-2*rise_time]))/(2*I_0) #account for the resistor difference 10kohm ch1 5kohm chan2

    #FP_modulation_sensitiviy=FP_shift/FP_modulation_amplitude#3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
    d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second  #NB could also be determined in terms of the min and max wavelength
    Δν_h=-3e10 #1/Δim*(-3e9)#(average_ν-ν_0)#<v>-vo/deltap/p  

    println("#NB set (average_ν-ν_0) to constant and using scaling factors.")

    #setting up the sampling timing
    T=data_points_per_period#period in sample numbers 
    ts_range=Int.(range(rise_time,step=1,stop=31250)) #31250 is the number of data points per half modulation cycle# important in the averaging functions. 
    averages=Int(length(ts_range)) #number of data points to average over in each modulation cycle
   
    #calculate the frequencyc/wavelength scale from FP shift.

    d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second
    FP_min_wavelength=c/ν_0-(FP_shift/2)
    λ=FP_min_wavelength .+ d_lambda_d_t .* range(0,step=0.001,length=sample_periods)
    ν=transpose([c]./[λ;])

    #set the pulse start values i.e triggers
    #indices=trigger_values(LD,T,sample_periods) #
    indices=triggers.-offset.+1
    ##get the 'direct' or reference data i.e. no modulation
    L_direct=matread(fp_direct_file)["data"]
    I0_direct=matread(ld_direct_file)["data"]
    I0_direct=mean(I0_direct)*resistor_value
    L_direct=L_direct/I0_direct #normalize the direct measurement. NB should use its own reference channel if possible. since current setting may not be exactly I_0
    L_direct=L_direct[offset:end]  #there is no reason for the offset rather than to have the data the same length as the modulated data. 
    L_direct_avg=LD_average_values(L_direct,indices,sample_periods,T) 
    L_direct_avg=L_direct_avg[1:sample_periods] 

    #calculate the sum and difference signals
    i_plus, i_minus = average_values(FP,indices,sample_periods,averages,T,ts_range,I_0) # No change needed

    # Calculate L and L_prime functions
    L = L0.(i_plus, i_minus, Δim) # Assuming Δim is available
    L_prime = L0_prime.(i_plus, i_minus, Δν_h, Δim) # Assuming Δν_h and Δim are available

    #divide the data into pulses 
    #these are offset by 1 to allow for subtracting the rise time from the first pulse.
    FP_array, LD_array = make_arrays(FP,LD,sample_periods-1,T,indices[2:end],rise_time)
            #old
            #L_direct_fit, L_prime_direct = lorentzian_fit(λ,L_direct_avg)
            #L_prime_direct = convert_to_dv(L_prime_direct,λ)

    #use a lorentzian fit of the direct data to create a reference derivative function L_direct_prime.
    #scale the calculated derivative to this value using a fit. The error signal is the difference in these two.
    scale_factor=0
    L_prime_direct=0
    try #for the stability files where the fit will fail.
        scale_factor,L_prime_direct=calculate_derivative_error(λ,L_prime,L_direct_avg)

        fwhm_method, fwhm_direct = get_fwhm(λ, L, L_direct_avg)
        println("FWHM method: $fwhm_method, FWHM direct: $fwhm_direct")
    catch
        println("If this is not a stability file,something is wrong.")
    end

   return λ, ν, L, L_prime, L_direct_avg, L_prime_direct,FP_array, LD_array,scale_factor,Δim
end
function load_and_return_raw_data(fp_file,ld_file,fp_direct_file,ld_direct_file)
    FP_data = matread(fp_file)
    FP = FP_data["data"]
    time = FP_data["time"]
    LD = matread(ld_file)["data"]
    FP_data = nothing  # Free memory

    # Use experiment parameters from exp_params dictionary
    rise_time=250 #the rise time of the intensity signal from the trigger point. total risetime = 2x risetime
    data_points_per_period=Int(0.001/sample_interval)
    T=data_points_per_period
    sample_periods=999
    ## If data triggered on FP, allign the LD pulses and adjust data accordingly. 
    LD=LD*resistor_value
    LD_mean=mean(LD) #find the mean value to set the 'trigger' point. 
    if LD_mean ==Inf
        error("infinite error")
    end
    triggers=Int.(trigger(LD,T,sample_periods))
    offset=triggers[1]
    LD=LD[offset:end] 
    FP=FP[offset:end]
    time=time[offset:end]
    #derived values
   return FP,time
end

function L0(i_p,i_m,Δim)
    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m,Δν_h,Δim)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function average_values(FP,indices,sample_periods,averages,T,ts_range,I_0)
    i_plus=zeros(sample_periods,averages)
    i_minus=zeros(sample_periods,averages)
    FP_array=zeros(T,sample_periods)
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

#treat the direct LD the same as the modulated data, except average for each period. 
function LD_average_values(L_direct,indices,sample_periods,T)
    #should fit this function instead of averaging
    LD_array=zeros(sample_periods)
    for i in 1:sample_periods
        LD_array[i]=mean(L_direct[indices[i]:Int(indices[i]+T-1)])
    end
   return LD_array
end

function dλ_to_dν(dI_dλ,x)
    c=2.99e8
    λ= x.*1e-9 
    dI_dν=zeros(length(dI_dλ))

    dλ_dν=(λ.^2)./c
    ν=c./(λ)
    #c/λ^2
    for i in 1:length(dλ_dν)
        dI_dν[i]=dI_dλ[i]*1e-9*dλ_dν[i] #1e-9 to convert the derivative in units of nm to units of meters. 
    end
    return dI_dν, ν
end
function convert_to_dv(L_prime_direct,λ)
    #convert to derivative to dI/dν from dI/d_lambda
    #direct_derivative, ν = dλ_to_dν(L_direct_fit,x) #fit params and lambda in nm 

    dI_dν=zeros(length(L_prime_direct))
    dλ_dν=(λ.^2)./c
    ν=c./(λ)
        #c/λ^2
    for i in 1:length(dλ_dν)
            dI_dν[i]=L_prime_direct[i]*1e9*dλ_dν[i] #1e-9 to convert the derivative in units of nm to units of meters. 
    end
    return dI_dν
end
#this really just prepares the data by formatting it the same way and dividing into each pulse.
function make_arrays_old(FP,LD,sample_periods,T,indices)
    FP_array=zeros(T,sample_periods)
    direct_array=zeros(T,sample_periods)
    LD_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        LD_array[:,i]=LD[Int(indices[i]):Int(indices[i]+T-1)] 
        FP_array[:,i]=FP[Int(indices[i]):Int(indices[i]+T-1)]./LD_array[:,i] #FP/LD
        #direct_array[:,i]=L_direct[indices[i]:Int(indices[i]+T-1)] use L direct_avg
    end
    return FP_array,LD_array
    
end
#starts the pulse at the initial rise, not the mid point, i.e normal trigger.Could cause problems on the first pulse, if triggered on the LD and not the function generator..
function make_arrays(FP,LD,sample_periods,T,indices,rise_time)
    FP_array=zeros(T,sample_periods)
    LD_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        LD_array[:,i]=LD[Int(indices[i]-rise_time):Int(indices[i]+T-1-rise_time)] 
        FP_array[:,i]=FP[Int(indices[i]-rise_time):Int(indices[i]+T-1-rise_time)]./LD_array[:,i] #FP/LD
        #direct_array[:,i]=L_direct[indices[i]:Int(indices[i]+T-1)] use L direct_avg
    end
    return FP_array,LD_array
    
end
function get_chirp_value(chirp, x,L_direct_avg,λ)
    # Check if x is within the bounds of the array 'L_direct_avg'
    if x < 1 || x > length(L_direct_avg)
        error("Starting point x is out of bounds in vector L_direct_avg")
    end
    chirp=rollmean(chirp,100)
    #should average first. 
    max_value,max_position=findmax(chirp)
    min_value,min_position=findmin(chirp)
    # Determine which one (max or min) occurs first in 'chirp'
    if max_position < min_position
        y, z = max_value, min_value
    else
        y, z = min_value, max_value
    end

    # Initialize indices for y and z in vector 'L_direct_avg'
    index_of_y = 1
    index_of_z = 1
    # Search backward for y in 'L_direct_avg'
    for i in x:-1:1
        if abs(L_direct_avg[i] - y) <= .0035 #calculated from spacing in L_direct_avg.
            index_of_y = i
            break
        end
    end

    # Search forward for z in 'L_direct_avg'
    for i in x:length(L_direct_avg)
        if abs(L_direct_avg[i] - z) <= .0035
            index_of_z = i
            break
        end
    end

    return abs(λ[index_of_y]-λ[index_of_z])
end

function get_chirp_estimate(FP_array, LD_array,L_direct_avg,λ, N)
    chirp=zeros(N)
    for i in 150:N-150
        chirp[i]=get_chirp_value(FP_array[:,i],i,L_direct_avg,λ)
    end
    return chirp
end

function lorentzian_fit(λ,L_direct_avg)
        ################################################################ lorentzian fit 
   # @. lor_model(x,p)= p[3]*(1/2*p[2])/(pi*(x-p[1])^2+(1/2*p[2])^2)#lorentzian model
    #@. lor_derivative_model(x,p)= p[3]*(-16(x-p[1])*p[2])/(pi*(4(x-p[1])^2+p[2]^2))^2#lorentzian model
    x=λ*1e9
    y=L_direct_avg
    p0=[1547,0.2,1.0] #initial guess
    lor_fit=curve_fit(lor_model,x,y,p0)

    L_direct_fit=lor_model(x,lor_fit.param)# 
    L_prime_direct=lor_derivative_model(x,lor_fit.param) #dI/dlambda nm
    return L_direct_fit, L_prime_direct
end



function get_fwhm(λ, L,L_direct_avg)
        ################################################################# 
    #calculate the fwhm
    ind1=findfirst(x->x>(maximum(L)/2),L)
    ind2=findfirst(x->x==maximum(L),L)
    fwhm_method=abs(2*(λ[ind1]-λ[ind2]))
    ind1=findfirst(x->x>(maximum(L_direct_avg)/2),L_direct_avg)
    ind2=findfirst(x->x==maximum(L_direct_avg),L_direct_avg)
    fwhm_direct=abs(2*(λ[ind1]-λ[ind2]))
    return fwhm_method, fwhm_direct
end

# function trigger_values(LD,T,sample_periods,search_width=1000)
#     LD_trigger=LD.-mean(LD)
#     #search_width=1000 #based on visually checking first and last LD pulse
#     step_size=Int(T/2-search_width/2)
#     indices=zeros(sample_periods*2)
#     search_start=1
#     for i in 1:length(indices) 
#         next_ind=findmin(abs.(LD_trigger[search_start:search_start+step_size]))[2]#find the minimum element in the expected range 
#         indices[i]=Int(search_start+next_ind)
#         search_start=search_start+next_ind+step_size
#     end
#     indices=Int.(indices[1:2:(2*sample_periods)])
#     return indices
# end

# function offset_trigger_values(LD,T,sample_periods,search_width=1000)
#     LD_trigger=LD.-mean(LD)
#     #search_width=1000 #based on visually checking first and last LD pulse
#     step_size=Int(T/2-search_width/2)
#     indices=zeros(sample_periods*2)
#     search_start=1

#     for i in 1:length(indices) 
#         next_ind=findmin(abs.(LD_trigger[search_start:search_start+step_size]))[2]#find the minimum element in the expected range 
#         indices[i]=Int(search_start+next_ind)
#         search_start=search_start+next_ind+step_size
#     end
#     #indices=Int.(indices[1:2:(2*sample_periods-1)]) removed.
#     return Int.(indices)

# end

function trigger(LD, T, sample_periods)
    # Offset LD by its mean
    LD_trigger = LD .- mean(LD) #center about zero

    # Step 1: Find the first rising trigger
    half_period = Int(floor(T / 2))
    
    # Find the minimum value in the first half-period
    min_value = minimum(abs.(LD_trigger[1:half_period]))
    min_index = findfirst(x -> abs.(x) == min_value, LD_trigger[1:half_period])
    # Check the slope around this trigger point using plus/minus 5 data points
    first_zero=0
    i=1 #loop counter
    rising_trigger=false
    while rising_trigger==false
        # Find the minimum value in the first half-period
        min_value = minimum(abs.(LD_trigger[1+(i-1)*half_period:half_period+(i-1)*half_period]))
        min_index = (i-1)*half_period+findfirst(x -> abs.(x) == min_value, LD_trigger[1+(i-1)*half_period:half_period+(i-1)*half_period])
        if min_index > 5 && min_index + 5 <= length(LD_trigger)
            if mean(LD_trigger[min_index-5:min_index]) < mean(LD_trigger[min_index:min_index+5])
                # It's a rising trigger point
                first_zero=min_index
                println("first zero =$first_zero")
                rising_trigger=true
            end
        else
            println("Trigger point is too close to the start or end")
        end
        i=i+1
        if i > 2*sample_periods
            error("while loop error")
        end
    end

    indices=zeros(sample_periods)
    indices[1]=Int(first_zero) 

    # Proceed with finding the minimum of the next periods
   for i in 2:sample_periods
             # Reference the last index to avoid drift
             period_start = indices[i-1]  + T - 500  # 1000 data points around one period after the last trigger
             period_end = indices[i-1]  + T + 500
             # Adjust if the range goes beyond the length of LD_trigger
             period_start = Int(max(1, period_start))
             period_end = Int(min(length(LD_trigger), period_end))

            period_min_value = minimum(abs.(LD_trigger[period_start:period_end]))
            period_min_index = findfirst(x -> abs.(x) == period_min_value, LD_trigger[period_start:period_end]) + period_start - 1

            indices[i]=Int(period_min_index)
    end
    #check spacing
    check=diff(indices)
    if maximum(check)>62570 || minimum(check) <62430
        error("the trigger points are not equally spaced, significant timing error greater than 70*16 ns")
    end

    return Int.(indices)
end


# function find_rising_trigger(LD, indices, window_size)
#     for index in indices
#         ##check the slope around the point
#         if index>window_size 
#             if LD[index-window_size]<LD[index+window_size]
#                 return index
#             end    
#         end
#     end
#     return -1
# end

#####plots etc
# function multi_plot()
#     # Initialize the plots
#     plot_L = plot(title = "L")
#     plot_L_prime = plot(title = "L_prime")
#     plot_L_direct_avg = plot(title = "L_direct_avg")
#     plot_L_prime_direct = plot(title = "L_prime_direct")

#     for (index, (λ, ν, L, L_prime, L_direct_avg, L_prime_direct)) in pairs(results)
#         # Add each series to the corresponding plot
#         plot!(plot_L, λ, L, label="L $index")
#         plot!(plot_L_prime, λ, L_prime, label="L_prime $index")

#         plot!(plot_L_direct_avg, λ, L_direct_avg, label="L_direct_avg $index")
#         plot!(plot_L_prime_direct, λ, L_prime_direct, label="L_prime_direct $index")
#     end
#     display(plot_L)
#     display(plot_L_prime)
#     display(plot_L_direct_avg)
#     display(plot_L_prime_direct)
# end

function multi_plot(results,scale_factor,modulation,save)
    # Initialize the plots
    plot_L = plot(xlabel="shift (pm)",ylabel="transmittance",dpi=600)
    plot_L_prime = plot(xlabel="shift (pm)",ylabel="derivative of transmittance (nm/Hz))",dpi=600)
    plot_L_direct_avg = plot(title = "L_direct_avg",xlabel="shift (pm)",ylabel="transmittance",dpi=600)
    plot_L_prime_direct = plot(title = "L_prime_direct",xlabel="shift (pm)",ylabel="derivative of transmittance (nm/Hz))",dpi=600)
    
    for (index, (λ,λ2,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array)) in pairs(results)
        mod=round(modulation[index],digits=2)
        # Add each series to the corresponding plot
        plot!(plot_L, λ*1e12, L, label="ΔI/I $mod")
        plot!(plot_L_prime, λ*1e12, L_prime*scale_factor[index], label="$mod")

        plot!(plot_L_direct_avg, λ2*1e12, L_direct_avg, label="L_direct_avg $index")
        plot!(plot_L_prime_direct, λ2*1e12, L_prime_direct, label="L_prime_direct $index")
    end
    display(plot_L)
    display(plot_L_prime)
    display(plot_L_direct_avg)
    display(plot_L_prime_direct)
    if save==true
        savefig(plot_L,"L.png")
        savefig(plot_L_prime,"L_prime.png")
    end
    #savefig(plot_L_direct_avg)
    #savefig(plot_L_prime_direct)
end

# function plot_max_L()
#     # Extracting maximum L values and corresponding indices
#     max_L_values = []
#     indices = []

#     for (index, (λ, ν, L, L_prime, L_direct_avg, L_prime_direct)) in pairs(results)
#         push!(max_L_values, maximum(L))
#         push!(indices, index)
#     end

#     # Sorting the values based on indices for a coherent plot
#     sorted_indices= sortperm(indices)

#     # Plotting
#     p=plot(indices[sorted_indices],max_L_values[sorted_indices], title="Maximum L vs Index", xlabel="Index", ylabel="Maximum L", label="Max L", marker=:circle)
#     display(p)
# end
function plot_max_L(results)
    # Extracting maximum L values and corresponding indices
    max_L_values = []
    indices = []

    for (index, ( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array)) in pairs(results)
        push!(max_L_values, maximum(L))
        push!(indices, index)
    end

    # Sorting the values based on indices for a coherent plot
    sorted_indices= sortperm(indices)

    # Plotting
    p=plot(indices[sorted_indices],max_L_values[sorted_indices], title="Maximum L vs Index", xlabel="Index", ylabel="Maximum L", label="Max L", marker=:circle)
    display(p)
end

# Objective function: Sum of squared differences
function objective_function(x, L_prime, L_prime_direct)
    return sum((x .* L_prime .- L_prime_direct).^2)
end

# Optimization algorithm to find the best scaling factor
function find_scaling_factor(L_prime, L_prime_direct, initial_x=1.0, delta_x=0.001, tolerance=1e-25, max_iter=100000)
    #The tolerance needs to be set correctly to be meaningful. Could really just get rid of it. 
    x = initial_x
    f_min = objective_function(x, L_prime, L_prime_direct)
    for i in 1:max_iter
        f_plus = objective_function(x + delta_x, L_prime, L_prime_direct)
        f_minus = objective_function(x - delta_x, L_prime, L_prime_direct)
        if f_plus < f_min
            x += delta_x
            f_min = f_plus
        elseif f_minus < f_min
            x -= delta_x
            f_min = f_minus
        else
            if abs(f_min) < tolerance
                println("Convergence reached with tolerance at iteration: ", i," preffer if no further improvements is shown instead?")
            else
                println("No further improvement at iteration with this delta_x", i)
            end
            break
        end
    end
    return x
end


function get_scaling_factors(results,modulation_values)
    #scale the derivative by matching to to the direct derivative. Which is itself a fit to the measured direct transmission.
    scaling_factors = []

    for index in 1:length(modulation_values)
        _, _, _, L_prime, _, L_prime_direct = results[modulation_values[index]]
        scaling_factor = find_scaling_factor(L_prime, L_prime_direct)
        push!(scaling_factors, scaling_factor)

    end
    
    return scaling_factors
end

# Function to find scaling factors for each modulation index
# function find_normalized_scaling_factors(results)
#     #scale the derivative by matching to to the direct derivative. Which is itself a fit to the measured direct transmission.
#     scaling_factors = []

#     for index in 1:length(modulation_values)
#         _, _, _, L_prime, _, L_prime_direct = results[modulation_values[index]]
#         scaling_factor = find_scaling_factor(L_prime, L_prime_direct)
#         push!(scaling_factors, scaling_factor)

#     end
#     scaling_factors=[1]./(scaling_factors./(scaling_factors[end]))
#     return scaling_factors
# end

function derivative_scaling_visual_check(results,scaling_factor)
    
    plot_L_prime_direct = plot(title = "L_prime_direct")
    i=1
    for (index, (λ, ν, L, L_prime, L_direct_avg, L_prime_direct)) in pairs(results)
        # Add each series to the corresponding plot
        plot!(plot_L_prime_direct, λ, scaling_factor[i]*L_prime, label="L_prime_direct")
        plot!(plot_L_prime_direct, λ, L_prime_direct, label="L_prime_direct $index")
        i=i+1
    end
   
    display(plot_L_prime_direct)
end

function save_results(results, filename)
    open(filename, "w") do file
        serialize(file, results)
    end
    println("Results saved to $filename")
end

function load_results(filename)
    open(filename, "r") do file
        return deserialize(file)
    end
end

# function simple_chirp_estimate(FP_array,L_direct_avg,λ, window) #window in nm, lambda in m
#     #calculate the chirp over a portion of the recorded signal defined by the window size and the peak position. 
#     #shouldn't be used on flat data, i.e stability measure
#     λ=λ*1e9
#     max_index=findfirst(x->x==maximum(L_direct_avg),L_direct_avg)
#     min_index=findfirst(x->abs(x-λ[1])>window,λ)
#     L_direct_avg
#     chirp=zeros(length(L_direct_avg))
#     FP=FP_array[:,min_index:max_index-min_index]
#     for i in 1:length(FP[1,:])
#         foo= exponential_fit(rollmean(FP[31500:end-250,i],200))
#         #FP[250:31000,i]#rollmean(FP_array[250:31000,i],200)
#         idx1=min_index+findfirst(x->x>=maximum(foo),L_direct_avg)
#         idx2=min_index+findfirst(x->x>=minimum(foo),L_direct_avg)
#         chirp[i]=abs(λ[idx1]-λ[idx2])
#     end
#    return chirp
# end

# function simple_chirp_estimate_stab(FP,L_direct_avg,λ) 
#     λ=λ*1e9
#     L_direct_avg
#     chirp=zeros(length(FP_array[1,:]))
#     for (i, slice) in enumerate(eachslice(FP, dims=2))
#         high=  exponential_fit(slice[250:31000])
#         low= exponential_fit(slice[31500:end-250])
#         #foo=(rollmean(FP[:,i],200))
#         #foo2=rollmean(FP[250:31000,i],200)
#         #FP[250:31000,i]#rollmean(FP[250:31000,i],200)
#         idx1=findfirst(x->x>=maximum(high),L_direct_avg)
#         idx2=findfirst(x->x>=minimum(low),L_direct_avg)
#         chirp[i]=abs(λ[idx1]-λ[idx2])
#     end
#    return chirp
# end


# function exponential_fit(FP)
    
#     ################################################################ lorentzian fit 
#     @. exp_model(x,p)=p[3]+p[1]*(1-exp(-x/p[2]))#lorentzian model
#     p0=[0.25,1/length(FP)] #initial guess
#     y=FP
#     x=1:1:length(FP)
#     chirp_fit=curve_fit(exp_model,x,y,p0)
#     chirp=exp_model(x,chirp_fit.param)# 
#     return chirp,chirp_fit
# end
# @. function lor_model(x,p)

    # p[3]*(1/2*p[2])/(pi*(x-p[1])^2+(1/2*p[2])^2)#lorentzian model
    # end
function lorentzian_fit_params(λ,L_direct_avg)
    x=λ*1e9
    y=L_direct_avg
    p0=[1547,0.2,1.0] #initial guess
    lor_fit=curve_fit(lor_model,x,y,p0)

return lor_fit.param
end

# function exponential_fit(FP,L_direct_avg,λ)
#     lor_params=lorentzian_fit_params(λ,L_direct_avg) #nb lambda converted to nm
#     x=λ*1e9 # to use same units in this function as the fit.
#     p_l=lor_params
#     @. lor_model(x)= p_l[3]*(1/2*p_l[2])/(pi*(x-p_l[1])^2+(1/2*p_l[2])^2)#lorentzian model
#     function modified_lor_model(x)
#         lambda0= 1546.647553073142
#         if x > lambda0
#             return p_l[3] * (1/2 * p_l[2]) / (pi * (x - p_l[1])^2 + (1/2 * p_l[2])^2)
#         else
#             return 0
#         end
#     end
    



#     @. exp_model(t,p)=modified_lor_model.(p[3]+p[1]*(1-exp(-t/p[2]))+p[4]t)#lorentzian model
#     B=.2 #pm/ma
#     τ=3000#microsecond
#     λ0=1546.5

#     p0=[B,τ,λ0,1]
    
#     #initial guess
#     t=range(start=0,step=1,length=length(FP)) #not actually time just sample number
#     chirp_fit=curve_fit(exp_model,t,FP,p0)
#     chirp=exp_model(t,chirp_fit.param)# 
#     return chirp,chirp_fit
# end

###########functions to account for drift in the peak positions of the measurements
function align_and_trim_vectors(vectors)
    function find_closest_to_zero(v)
        third_len = length(v) ÷ 3
        center_start = Int(round(third_len + 1))
        center_end = Int(round(center_start + third_len))
        center_segment = v[center_start:center_end]
        _, min_index = findmin(abs.(center_segment))
        return Int(round(center_start + min_index - 1))
    end

    zero_indices = [find_closest_to_zero(v) for v in vectors]
    distances_to_ends = [(length(v) - zero_idx, zero_idx - 1) for (v, zero_idx) in zip(vectors, zero_indices)]
    min_distance_to_end = minimum([min(d[1], d[2]) for d in distances_to_ends]) #

    n_new = 2 * min_distance_to_end + 1  # New length based on the shortest distance to end
    aligned_vectors = []
    for (v, zero_idx) in zip(vectors, zero_indices)
        start_idx = max(1, Int(round(zero_idx - min_distance_to_end)))
        end_idx = min(length(v), Int(round(zero_idx + min_distance_to_end)))
        println(start_idx,end_idx)
        adjusted_v = v[start_idx:end_idx]
        push!(aligned_vectors, adjusted_v)
    end

    return aligned_vectors
end

function center_pairs_in_dict(results,keys)
    #chatgpt assisted
    keys_of_interest = keys
    new_results = deepcopy(results) # Create a copy of the results dictionary
    
    #choose a set of vectors corresponding to one dataset, eg. L,L_prime
    for (prime_index, associated_index) in [(4, 3), (6, 5)]
        prime_vectors = [new_results[key][prime_index] for key in keys_of_interest if haskey(new_results, key)]
        associated_vectors = [new_results[key][associated_index] for key in keys_of_interest if haskey(new_results, key)]

        centered_prime_vectors = align_and_trim_vectors(prime_vectors)
        
        #update all the datasets in results, with the new vectors.
        for (i, key) in enumerate(keys_of_interest)
            if haskey(new_results, key)
                current_entry = new_results[key]
                start_idx, end_idx = get_centering_indices(centered_prime_vectors[i], prime_vectors[i])
                centered_associated_vector = associated_vectors[i][start_idx:end_idx]

                updated_entry = (
                                 prime_index == 4 ? current_entry[1][1:length(centered_prime_vectors[i])] : current_entry[1], #these lines update lambda and nu to be lambda1 and lambda2 corresponding to the lengths of the direct and regular datasets after centering.
                                 prime_index == 6 ? results[key][1][1:length(centered_prime_vectors[i])]  : current_entry[2],
                                 prime_index == 4 ? centered_associated_vector : current_entry[3],
                                 prime_index == 4 ? centered_prime_vectors[i] : current_entry[4],
                                 prime_index == 6 ? centered_associated_vector : current_entry[5],
                                 prime_index == 6 ? centered_prime_vectors[i] : current_entry[6],
                                 current_entry[7], current_entry[8])
                new_results[key] = updated_entry
            end
        end
    end

    return new_results
end

function center_lambda(results,keys)
    center_index_1=findmin(abs.(results[keys[end]][4]))[2]
    center_index_2=findmin(abs.(results[keys[end]][6]))[2]
    λ=results[keys[1]][1]
    λ1=λ.-λ[center_index_1]
    λ2=results[keys[1]][2].-λ[center_index_2]

    new_results=deepcopy(results)
    for i in keys
        updated_entry = (
            λ1,
            λ2,
            new_results[i][3],
            new_results[i][4],
            new_results[i][5],
            new_results[i][6],
            new_results[i][7],
            new_results[i][8])

        new_results[i]=updated_entry
    end
    return new_results
end


function get_centering_indices(centered_vector, original_vector)
    # Assuming centered_vector is a segment of original_vector
    start_idx = findfirst(isequal(centered_vector[1]), original_vector)
    end_idx = findlast(isequal(centered_vector[end]), original_vector)
    return start_idx, end_idx
end

function get_intensity_modulation(results,keys)
    modulation_dict=Dict()
    modulation_vect=zeros(length(keys))
    for (i, key) in enumerate(keys)
        if haskey(results, key)
            I=results[key][8]
            I0=mean(mean(I,dims=1))
            Ih=mean(mean(I[32000:end-750,:],dims=1)) # this will depend on the trigger point, but these are large buffers.
            Il=mean(mean(I[750:30500,:],dims=1))
            mod=abs(Ih-Il)/(2*I0)
            
            if abs(abs(I0-Ih)/I0-mod)>0.01
               println("this should never happen")
               break
            end
            modulation_dict[key] = mod
            modulation_vect[i] =  mod    
        end
    end
    return modulation_dict,modulation_vect
end

# Example usage:
# Assuming results is your dictionary and align_and_trim_vectors is defined as before

function multi_chirp_estimate(start,stop,results,modulation_values)
    chirp=Dict()
    for i in modulation_values
        λ,_,_,_,L_direct_avg,_,FP_array, _=results[i] 
        chirp(i)=simple_chirp_estimate_stab(FP_array[:,start:stop],L_direct_avg,λ)
    end
    return chirp
end

####finding λ(t) without any fitting. 
# function direct_chirp_calculation(λ, L_direct_avg,FP)
#     center=findmax(L_direct_avg)[2] 
#     T1=L_direct_avg[1:center]  # can improve the resolution by using the lorentzian fit and an arbitrary lambda scale.
#     T2=L_direct_avg[center:end] 
#     chirp1=zeros(length(FP))
#     chirp2=zeros(length(FP))
#     for i in 1:length(chirp1)
#        idx =findmin(abs.(FP[i].-T1))[2] ## 
#        chirp1[i]=λ[Int(idx)]  
#     end
#     for i in 1:length(chirp1)
#         idx =findmin(abs.(FP[i].-T2))[2] 
#         chirp2[i]=λ[Int(idx)] 
#     end
       
#     return chirp1,chirp2 
# end
function direct_chirp_calculation(λ, L_direct_avg,FP)
    fit_params = lorentzian_fit_params(λ, L_direct_avg) #the function converts lambda to nm
    λ_0=fit_params[1]*1e-9 
    λ=range(0,step=1e-13,stop=0.2e-9).+λ_0
    L=lor_model(λ*1e9,fit_params)
    chirp=zeros(length(FP))
    for i in 1:length(chirp)
       idx =findmin(abs.(FP[i].-L))[2] ## 
       chirp[i]=λ[Int(idx)]  
    end
    return -chirp,0 #_ to be compatible with previous function
end

function calculate_derivative_error(λ, L_prime, L_direct_avg)

    #nb the fit is done in nm no m

    fit_params = lorentzian_fit_params(λ, L_direct_avg) #the function converts lambda to nm

    # Estimating zero crossing
    idx1 = floor(Int, length(L_prime) / 3)
    idx2 = floor(Int, 2 * length(L_prime) / 3)
    min_val, min_idx = findmin(abs.(L_prime[idx1:idx2]))
    zero_crossing_idx = idx1 + min_idx - 1  # Adjusting index

    # Linear interpolation around the zero crossing
    n = 9  # Number of points to use for interpolation
    λ_interpolated = λ[max(1, zero_crossing_idx - n):min(length(λ), zero_crossing_idx + n)]
    L_prime_interpolated = L_prime[max(1, zero_crossing_idx - n):min(length(L_prime), zero_crossing_idx + n)]
    
    # Linear interpolation function
    interpolation_function = LinearInterpolation(λ_interpolated, L_prime_interpolated)

    # Finding zero of the interpolated function within tolerance T
    T = minimum(abs.(L_prime))/1000 # Set your tolerance
    λ_zero = find_zero(x -> interpolation_function(x), (λ_interpolated[1], λ_interpolated[end]), Bisection(), atol = T)
    # Updating λ in fit parameters
    
    fit_params[1] = λ_zero*1e9 #m*nm/m
    
    lorentzian_fit= lor_derivative_model(λ*1e9, fit_params) #1e9 nm
    lorentzian_fit=lorentzian_fit*1e9*λ_zero^2/2.99e8 #1/nm *nm/m *m^2/(m/s)=s = 1/hz =dI/dv ok.
    scale_factor=find_scaling_factor(L_prime,lorentzian_fit)
    #error = L_prime*scale_factor .- (lorentzian_fit)
    @infiltrate
    return scale_factor,lorentzian_fit 
end
function calculate_derivative_slope(λ, L_direct_avg)
    fit_params = lorentzian_fit_params(λ, L_direct_avg) #the function converts lambda to nm
    second_derivative= lor_second_deriv_model(fit_params[1], fit_params) #1e9 nm
    λ_zero=fit_params[1]*1e-9 #nm*m/nm

    second_derivative=second_derivative*1e9^2#*λ_zero^2/2.99e8 #1/nm² *nm/m *m^2/(m/s)=s = 1/hz =dI/dv ok.
       return second_derivative
end
function get_resonator_length(lambda,L_direct_avg,num_oscillations,min_value,max_value,rise)
    lambda1=0
    lambda2=0
    if rise==false
        lambda1=findfirst(x->x<min_value,L_direct_avg[500:end])
        lambda2=findfirst(x->x<max_value,L_direct_avg[500:end])
        
    elseif rise==true
        lambda1=findfirst(x->x>min_value,L_direct_avg[:])
        lambda2=findfirst(x->x>max_value,L_direct_avg[:])
    end
    delta_lambda=abs(lambda[lambda1]-lambda[lambda2])/num_oscillations
    delta_f=3e8*(delta_lambda)/(1547e-9^2)
    n=1.46
    resonator_length=3e8/(n*2*delta_f)
end

# function calculate_derivative_error(λ,L_prime,L_direct_avg)
#     @. lor_derivative_model(x,p)= p[3]*(-16(x-p[1])*p[2])/(pi*4(x-p[1])^2+p[2]^2)^2#lorentzian model
#     fit_params=lorentzian_fit_params(λ,L_direct_avg)

#     #in the function f(\lambda)=L_prime, find f(\lambda) = zero. Do this by interpolation the functions.

#     fit_params[1]= #new λ.
#     lorentzian_fit=lor_derivative_model(λ,fit_params)
#     error= L_prime.-lorentzian_fit
#     return λ,L_prime,L_prime_direct 
# end

#p1=x0,p[2]=gamma p[3] is scaling
@. lor_model(x,p)=p[3]/pi*  0.5*p[2]/((x-p[1])^2+(0.5*p[2])^2)
@. lor_derivative_model(x,p)= p[3]*(-16(x-p[1])*p[2])/(pi*(4(x-p[1])^2+p[2]^2)^2)#lorentzian model
@. lor_second_deriv_model(x,p)=p[3]*16*p[2]*(12*(x-p[1])^2-p[2]^2)/(pi*(4(x-p[1])^2+p[2]^2)^3)
