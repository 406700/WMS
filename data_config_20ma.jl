
dir="/home/michael/Documents/Experimental_Data/WMS_final_paper/20240104/"
fp_file =  dir*"20mv_ch1.mat"
ld_file = dir* "20mv_ch2.mat"#"path/to/first/ld/file.mat"

#FP 
#FP 
FP_frequency=0.5 #hertz
FP_shift=3.63*1e-9 #nm*volt *1e-9m/nm
FP_modulation_amplitude=0.3 #volts
FP_shift=(FP_modulation_amplitude)*FP_shift

#LD 
c=2.99e8
ν_0=c/1546.637e-9


#direct files
#f_direct=readdir(dir)[3:4]
fp_direct_file=dir*"20mvnm_ch1.mat"
ld_direct_file=dir*"20mvnm_ch2.mat"


resistor_value=10/5 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)

