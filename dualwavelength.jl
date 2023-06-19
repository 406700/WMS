##NB possibly a slight drift in the signal over time, even though the trigger was set on the LD ref channel

using DelimitedFiles
using CSV
using DataFrames
using Plots
using Statistics
using MAT
plotlyjs()

#cd("/home/m/OneDrive/Experimental_Data/20230609")


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

function average_values_referenced(FP)
    #multiply the FP by the reference noise fluctuations (i.e scale the values)
    i_plus=zeros(sample_periods,averages)
    i_minus=zeros(sample_periods,averages)
    FP_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        FP_array[:,i]=FP[indices[i]:Int(indices[i]+T-1)]
    end
    #find the LD fluctuations
    LD_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        LD_array[:,i]=LD[indices[i]:Int(indices[i]+T-1)]
    end
    LD_mean=mean(LD_array[rise_time:end,:],dims=1)
    LD_noise=(LD_array)./LD_mean #find a scaling factor for I0
    LD_array=nothing #clear memory
    #return LD_noise
    for i in 1:averages
        ts=ts_range[i]
        i_plus[:,i]=(FP_array[ts,:]./LD_noise[ts,:]).+(FP_array[Int(ts+T/2),:]./LD_noise[Int(ts+T/2),:])
        i_minus[:,i]=(FP_array[ts,:]./LD_noise[ts,:]).-(FP_array[Int(ts+T/2),:]./LD_noise[Int(ts+T/2),:])
    end
    i_plus=mean(i_plus,dims=2)
    i_minus=mean(i_minus,dims=2)
    return i_plus/I_0,i_minus/I_0
end

####################################################CSV files
# cd("/home/m/OneDrive/Experimental_Data/20230605chirp_stability")
# data=empty
# f=readdir()[end]
# open(`head -n1000000 $f`) do io
#     global data=readdlm(io,',',skipstart=12)
# end
# header=empty
# open(`head -n12 $f`) do io
#     global header=readdlm(io,',')
#  end
# time=data[:,1].-data[1,1] #normalize to the trigger point
# FP=data[:,2] #Fabry perot transmission
# LD=data[:,3]*2 #Laser Diode reference times the resistor ratio. 2= R1/R2

############################################################# matlab files (Seems way faster, but maybe just from splitting up files for channels)

#cd("/home/m/OneDrive/Experimental_Data/20230609")
cd("/home/m/OneDrive/Experimental_Data/20230612_stability/")
f=readdir()#[3:4]#[3:4]
f=f[[10,12]]
FP_data=matread(f[1])
FP=FP_data["data"]
time=FP_data["time"]
LD=matread(f[2])["data"]
FP_data=nothing #free memory?

##if trigger not set to
resistor_value=2 #Rch1/Rch2
offset=1#54230#findfirst(x->x>0,time)
LD=LD[offset:end]*resistor_value
FP=FP[offset:end]
time=time[offset:end]
########################################################### parameters
#NB resistor values
sample_interval=1.6e-8 #seconds
data_points_per_period=Int(0.001/sample_interval)

rise_time=250 #risetime of intensity modulation. 250?
I_0=mean(LD[rise_time:data_points_per_period])
Δim= (mean(LD[rise_time:31250])-mean(LD[rise_time+31250:data_points_per_period]))/(2*I_0) #account for the resistor difference 10kohm ch1 5kohm chan2
FP_frequency=1 #hertz
measured_delta_lambda=7.2e-10 #my measure 720pm 820 gerogina
FP_modulation_amplitude=0.2 #volts
FP_modulation_sensitiviy=measured_delta_lambda/FP_modulation_amplitude#3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
FP_shift=FP_modulation_amplitude*FP_modulation_sensitiviy
d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second
Δν_h=1/Δim*(3e8/1546.920e-9-3e8/1546.896e-9)#<v>-vo/deltap/p  #NB delta im should be the one used in average measurement, but also this shouldn't change.

#getting the time
T=data_points_per_period#period in sample numbers 
ts=rise_time #defines the number of samples before the value stabilizes
ts_range=Int.(range(rise_time,step=1,stop=31250)) #31250 is the number of data points per half modulation cycle
sample_periods=999 #number of modulation cycles to average over
averages=Int(length(ts_range)) #number of data points to average over in each modulation cycle

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
################################################################ run code
i_plus,i_minus=average_values(FP)
i_plus_ref,i_minus_ref=average_values_referenced(FP)


L=L0.(i_plus,i_minus)
L_prime=L0_prime.(i_plus,i_minus)

#L_referenced=L0.(i_plus,i_minus)
#L_prime_referenced=L0_prime.(i_plus_ref,i_minus_ref)

###################################################################### averaging to I_0 for each pulse
    ##deleted something here
    #calculated the average I_0 and \delta Im for each pulse
    # I_0_t=zeros(sample_periods)
    # Δim_t=zeros(sample_periods)
    # for i in 1:sample_periods
    #     I_0_t[i]=mean(LD[indices[i]+rise_time:data_points_per_period+indices[i]])
    #     Δim_t[i]= (mean(LD[indices[i]+rise_time:indices[i]+31250])-mean(LD[indices[i]+rise_time+31250:indices[i]+data_points_per_period]))/(2*I_0_t[i])
    # end
    # L_referenced=L0_referenced.(i_plus,i_minus,I_0_t,Δim_t)
    # L_prime_referenced=L0_prime_referenced.(i_plus,i_minus,I_0_t,Δim_t)

    #pass deltaim to the L0 and prime functions. use ./ to divide I_plus_minus by I_0_t. 

################################################################  plots
#samples=range(1,step=1000,stop=28000)
plot(L,title="L")
plot(L_prime,title="L_prime")
#plot(L[:,samples],label="L") #for non averaged L
#plot(L_prime[:,samples],label="L_prime") #for non averaged L
#convert x axis from time to wavelength
t_to_lambda=range(0,step=.001,stop=0.001*sample_periods-.001)*d_lambda_d_t*1e9
plot(t_to_lambda,L)
################################################################# 
#FWHM 
ind1=findfirst(x->x>(maximum(L)/2),L)
ind2=findfirst(x->x==maximum(L),L)
fwhm=abs(2*(t_to_lambda[ind1]-t_to_lambda[ind2]))

plot(L)