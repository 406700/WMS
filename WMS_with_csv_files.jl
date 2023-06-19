#need to be careful of the trigger channel since the LD might not line up nicely with the trigger channel.
#Then an offset needs to be used to reallign the data. 
#Also potential for different resistor values in the setup.

using DelimitedFiles
using CSV
using DataFrames
using Plots
using Statistics
#using MAT
plotlyjs()

function L0(i_p,i_m)
    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function average_values(FP)
    i_plus=zeros(sample_periods,averages)
    i_minus=zeros(sample_periods,averages)
    for i in 1:averages
        ts=ts_range[i]
        sample_plus=range(ts,step=T,length=sample_periods)
        sample_minus=sample_plus.+Int(T/2)
        i_plus[:,i]=FP[sample_plus].+FP[sample_minus]
        i_minus[:,i]=FP[sample_plus].-FP[sample_minus]
    end
    i_plus=mean(i_plus,dims=2)
    i_minus=mean(i_minus,dims=2)
    return i_plus/I_0,i_minus/I_0
end

##############################################################CSV files
cd("/home/m/OneDrive/Experimental_Data/20230605chirp_stability")
data=empty
f=readdir()[end]
open(`head -n1000000 $f`) do io
    global data=readdlm(io,',',skipstart=12)
end
header=empty
open(`head -n12 $f`) do io
    global header=readdlm(io,',')
 end
time=data[:,1].-data[1,1] #normalize to the trigger point
FP=data[:,2] #Fabry perot transmission
LD=data[:,3]*2 #Laser Diode reference times the resistor ratio. 2= R1/R2

############################################################# matlab files 

#cd("/home/m/OneDrive/Experimental_Data/20230609")
# f=readdir()#[3:4]
# FP_data=matread(f[1])
# FP=FP_data["data"]
# time=FP_data["time"]
# LD=matread(f[2])["data"]
# FP_data=nothing #free memory?

# ##if trigger not set to
# resistor_value=2 #Rch1/Rch2
# offset=54230#findfirst(x->x>0,time)
# LD=LD[offset:end]*resistor_value
# FP=FP[offset:end]
# time=time[offset:end]
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
FP_modulation_sensitiviy=3.5e-9#meters/volt #measured calibration factor, not needed if measured_delta_lambda is available.
FP_shift=FP_modulation_amplitude*FP_modulation_sensitiviy
d_lambda_d_t=FP_shift/(1/(2*FP_frequency))#*1e12 #m #cycles/second
Δν_h=1/Δim*(3e8/1546.9205e-9-3e8/1546.896e-9)#<v>-vo/deltap/p  #NB delta im should be the one used in average measurement, but also this shouldn't change.

#getting the time
T=Int((1/1000)/1.6e-8) #period in sample numbers 
ts=rise_time #defines the number of samples before the value stabilizes
ts_range=Int.(range(rise_time,step=1,stop=31250)) #31250 is the number of data points per half modulation cycle
sample_periods=10 #number of modulation cycles to average over
averages=Int(length(ts_range)) #number of data points to average over in each modulation cycle

################################################################ run code
i_plus,i_minus=average_values(FP)

    #writedlm("iplus",i_plus)
    #writedlm("i_minus",i_minus)

    #i_plus=readdlm("iplus")
    #i_minus=readdlm("i_minus")

L=L0.(i_plus,i_minus)
L_prime=L0_prime.(i_plus,i_minus)

################################################################  plots
samples=range(1,step=1000,stop=28000)
plot(L,title="L")
plot(L_prime,title="L_prime")
plot(L[:,samples],label="L") #for non averaged L
plot(L_prime[:,samples],label="L_prime") #for non averaged L
#convert x axis from time to wavelength
t_to_lambda=range(0,step=.001,stop=0.001*sample_periods-.001)*d_lambda_d_t*1e9
plot(t_to_lambda,L)
################################################################# 
#FWHM 
#ind1=findfirst(x->x>(maximum(L)/2),L)
#ind2=findfirst(x->x==maximum(L),L)
#fwhm=abs(2*(t_to_lambda[ind1]-t_to_lambda[ind2]))
