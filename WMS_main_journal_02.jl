
using DelimitedFiles
using CSV
using DataFrames
using Plots
using Statistics
using MAT
using LsqFit
#using RollingFunctions
plotlyjs()

cd("/home/m/OneDrive/Experimental_Data/WMS_final_paper/")
cd("")

function L0(i_p,i_m)

    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function average_values(FP)
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
function LD_average_values(LD)
    LD_array=zeros(sample_periods)
    for i in 1:sample_periods
        LD_array[i]=mean(LD[indices[i]:Int(indices[i]+T-1)])
    end
   return LD_array
end

function dλ_to_dν(dI_dλ,x)
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

function chirp_estimate()
    #chirp_h = zeros(sample_periods)
    #chirp_l = zeros(sample_periods)
    FP_array=zeros(T,sample_periods)
    direct_array=zeros(T,sample_periods)
    LD_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        FP_array[:,i]=FP[indices[i]:Int(indices[i]+T-1)]./LD[indices[i]:Int(indices[i]+T-1)]
        LD_array[:,i]=LD[indices[i]:Int(indices[i]+T-1)]
        direct_array[:,i]=L_direct[indices[i]:Int(indices[i]+T-1)]
    end
    return FP_array,direct_array,LD_array
    
end

function get_chirp_value(chirp, x)
    # Check if x is within the bounds of the array 'L_direct_avg'
    if x < 1 || x > length(L_direct_avg)
        error("Starting point x is out of bounds in vector L_direct_avg")
    end

    # Find max and min values and their positions in vector 'chirp'
     max_value, max_position = findmax(chirp)
     min_value, min_position = findmin(chirp)

    # Determine which one (max or min) occurs first in 'chirp'
    if max_position < min_position
        y, z = max_value, min_value
    else
        y, z = min_value, max_value
    end

    #y=mean(chirp[30000:31000])
    #z=mean(chirp[61250:62250])
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

############################################################# matlab files (Seems way faster, but maybe just from splitting up files for channels) The

dir="/home/m/OneDrive/Experimental_Data/WMS_final_paper/"
fp_file =  dir*"03_ch1.mat"
ld_file = dir* "03_ch2.mat"#"path/to/first/ld/file.mat"

FP_data=matread(fp_file)
FP=FP_data["data"]
time=FP_data["time"]
LD=matread(ld_file)["data"]
FP_data=nothing #free memory?

#experiment parameters:

#can check if FP is > than LD, than must have had a higher resistance.
resistor_value=10/5 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)

#FP 
FP_frequency=0.5 #hertz
FP_min_wavelength=1546.286e-9-0.175e-9
FP_max_wavelength=1547.078e-9+0.175e-9
FP_shift= (FP_max_wavelength-FP_min_wavelength) #1546   my measure 720pm 820 gerogina.
FP_modulation_amplitude=0.3 #volts

#LD 
c=2.99e8
average_ν=c/1546.920e-9
ν_0=c/1546.896e-9
## If data triggered on FP, allign the LD pulses and adjust data accordingly. 
LD=LD*resistor_value
LD_mean=mean(LD) #find the mean value to set the 'trigger' point. 
offset=findfirst(x->x>=LD_mean,LD)[1] #middle of the first rise i.e trigger on positive pulse

LD=LD[offset:end] #adjust the voltage according to the resistor
FP=FP[offset:end]
time=time[offset:end]

########################################################### parameters


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
sample_periods=995 #number of modulation cycles to average over
averages=Int(length(ts_range)) #number of data points to average over in each modulation cycle
#λ=range(0,step=.001,stop=0.001*sample_periods-.001)*d_lambda_d_t .+FP_min_wavelength#0:.001:998 *dλ/dt *1e9 nm/m(NB nm/second)
λ=FP_min_wavelength .+ d_lambda_d_t .* range(0,step=0.001,length=995) #NB 499?
ν=transpose(c/λ)
#dI_to_dI/dν=1/d_lambda_d_t #dI/dt*1/d_lambda_d_t

#############################################################################
#create the timing and avoid drift in the LD (cant count on an exact number of samples per period)
LD_trigger=LD.-mean(LD)
search_width=1000 #based on visually checking first and last LD pulse
step_size=Int(T/2-search_width/2)
indices=zeros(sample_periods*2)
search_start=1
for i in 1:length(indices) 
    next_ind=findmin(abs.(LD_trigger[search_start:search_start+step_size]))[2]#find the minimum element in the expected range 
    indices[i]=Int(search_start+next_ind)
    global search_start=search_start+next_ind+step_size
end
indices=Int.(indices[1:2:(2*sample_periods-1)])

#unmodulated 
#the unmodulated data
fp_file =  dir*"03nm_ch1.mat"
ld_file = dir* "03nm_ch2.mat"#"path/to/first/ld/file.mat"

f=readdir()[3:4]
L_direct=matread(fp_file)["data"]
I_direct=matread(ld_file)["data"]*resistor_value

L_direct=L_direct/mean(I_direct) #normalize the direct measurement. NB should use its own reference channel if possible. since current setting may not be exactly I_0
#LD_direct_mean=matread(f[2])
L_direct=L_direct[offset:end]
L_direct_avg=LD_average_values(L_direct)
L_direct_avg=L_direct_avg[1:995] 

################################################################ lorentzian fit 
@. lor_model(x,p)= p[3]*(1/2*p[2])/(pi*(x-p[1])^2+(1/2*p[2])^2)#lorentzian model
@. lor_derivative_model(x,p)= p[3]*(-16(x-p[1])*p[2])/(pi*4(x-p[1])^2+p[2]^2)^2#lorentzian model
x=λ*1e9
y=L_direct_avg
p0=[1547,0.2,1.0] #initial guess
lor_fit=curve_fit(lor_model,x,y,p0)

L_direct_fit=lor_model(x,lor_fit.param)# 
L_prime_direct=lor_derivative_model(x,lor_fit.param) #dI/dlambda nm


################################################################ run code
i_plus,i_minus=average_values(FP)
L=L0.(i_plus,i_minus)
L_prime=L0_prime.(i_plus,i_minus)

L=L[1:995]
L_prime=L_prime[1:995]

#convert to derivative to dI/dν from dI/d_lambda
#direct_derivative, ν = dλ_to_dν(L_direct_fit,x) #fit params and lambda in nm 

dI_dν=zeros(length(L_prime_direct))
dλ_dν=(λ.^2)./c
ν=c./(λ)
    #c/λ^2
for i in 1:length(dλ_dν)
        dI_dν[i]=L_prime_direct[i]*1e9*dλ_dν[i] #1e-9 to convert the derivative in units of nm to units of meters. 
end

################################################################# 
#calculate the fwhm
ind1=findfirst(x->x>(maximum(L)/2),L)
ind2=findfirst(x->x==maximum(L),L)
fwhm=abs(2*(λ[ind1]-λ[ind2]))
ind1=findfirst(x->x>(maximum(L_direct_avg)/2),L_direct_avg)
ind2=findfirst(x->x==maximum(L_direct_avg),L_direct_avg)
fwhm=abs(2*(λ[ind1]-λ[ind2]))

##plotting comparison direct and indirect

plot(λ*1e9,L,label="L")
plot!(λ*1e9,L_direct_avg,label="L direct")
plot(x,dI_dν)
plot!(x,L_prime)


##plots for publication.
if false
    ##plotting the insert of the raw transmission data
    mid=findfirst(x->x>0.5,time)
    FPs=FP[mid:end]
    LDs=LD[mid:end]

    times=time[mid:end].-time[mid]
    max_index=findfirst(x->x>0.28,times)
    timeshift=findfirst(x->x>0.0999650,times) #tuned to start on the beggining of a cycle

    #plot the insert on left side. 
    tuning=
    start_index_l=max_index-timeshift
    stop_index_l=start_index_l+data_points_per_period*3
    plot(times[start_index_l:stop_index_l],FPs[start_index_l:stop_index_l])

    #plot the insert on right side. 
    #start_index_r=max_index+timeshift
    start_index_r=findfirst(x->x>.3800368,times) #tuning
    stop_index_r=start_index_r+data_points_per_period*3
    plot(times[start_index_r:stop_index_r],FPs[start_index_r:stop_index_r])

    #  matwrite("left_insert_30ma.mat", Dict(
    #         "time" => times[start_index_l:stop_index_l],
    #         "FP_voltage" => FPs[start_index_l:stop_index_l]
    #   ))

    #   matwrite("right_insert_30ma.mat", Dict(
    #     "time" => times[start_index_r:stop_index_r],
    #     "FP_voltage" => FPs[start_index_r:stop_index_r]
    # ))

    # #matwrite("left_insert_LD.mat", Dict(
    #         "time" => times[start_index_l:stop_index_l],
    #         "LD_voltage" => LDs[start_index_l:stop_index_l]
    #   ))

    #  # matwrite("right_insert_LD.mat", Dict(
    #     "time" => times[start_index_r:stop_index_r],
    #     "LD_voltage" => LDs[start_index_r:stop_index_r]
    # ))
    # #matwrite("FP.mat", Dict(
    #     "time" => times[:],
    #     "FP_voltage" => FPs[1:50:end]
    # ))
    # #matwrite("FP.mat", Dict(
    #     "time" => times[1:50:end],
    #     "FP_voltage" => FPs[1:50:end]
    # ))
    plot(FPs[1:10:end])

    #matwrite("L_30ma.mat", Dict(
    #    "wavelength_shift" => [t;f],
    #   "L" => L, 
    #  "L_prime"=> L_prime,
    #))
end


###chirp estimate 

FP_array,direct_array,LD_array=chirp_estimate()
chirp=zeros(995)
for i in 1:995
    chirp[i]=get_chirp_value(FP_array[:,i],i)
end
plot(chirp*1e12)
plot(λ ,L)
plot!(λ,L_direct_avg)
plot(L_prime)
plot!(dI_dν)