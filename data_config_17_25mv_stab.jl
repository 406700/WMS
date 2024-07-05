
#dir="/home/michael/Documents/Experimental_Data/WMS_final_paper/20240112_17ma_stab/"
dir="/home/michael/Documents/Experimental_Data/WMS_final_paper/20240114_stability/17/"

fp_file =  dir*"25_ch1.mat"
ld_file = dir* "25_ch2.mat"#"path/to/first/ld/file.mat"

#FP 
FP_frequency=0.5 #hertz
FP_shift=3.63*1e-9 #nm*volt *1e-9m/nm
FP_modulation_amplitude=0.3 #volts
FP_shift=(FP_modulation_amplitude)*FP_shift

#LD 
c=2.99e8
ν_0=c/1546.546e-9

#direct files

# #f_direct=readdir(dir)[3:4]
# fp_direct_file=dir*"25ref_ch1(0).mat"
# ld_direct_file=dir*"25ref_ch2(0).mat" #NB data file corrupted
fp_direct_file=dir*"ref_ch1.mat"    
ld_direct_file=dir*"ref_ch2.mat" 

resistor_value=10/10 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)

