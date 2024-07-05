
dir="/home/michael/Documents/Experimental_Data/WMS_final_paper/20240126/17/"

fp_file =  dir*"10_ch1.mat"
ld_file = dir* "10_ch2.mat"

#FP 
FP_frequency=0.5 #hertz
FP_shift=3.63*1e-9 #nm*volt *1e-9m/nm
FP_modulation_amplitude=0.3 #volts
FP_shift=(FP_modulation_amplitude)*FP_shift

#LD 
c=2.99e8
ν_0=c/1546.546e-9

#direct files
fp_direct_file=dir*"ref_ch1.mat"    
ld_direct_file=dir*"ref_ch2.mat"    

resistor_value=10/10 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)

