
dir="/home/m/OneDrive/Experimental_Data/WMS_final_paper/20231226/"
fp_file =  dir*"15mv_ch1.mat"
ld_file = dir* "15mv_ch2.mat"#"path/to/first/ld/file.mat"

#FP 
FP_frequency=0.5 #hertz
M=999
FP_min_wavelength=1546.286e-9-0.175e-9 #NB
FP_max_wavelength=1547.078e-9+0.175e-9 #NB
FP_shift= (FP_max_wavelength-FP_min_wavelength)*3/2 #1546   my measure 720pm 820 gerogina.
FP_modulation_amplitude=0.2 #volts

#LD 
c=2.99e8
average_ν=c/1546.920e-9
ν_0=c/1546.896e-9


#direct files
#f_direct=readdir(dir)[3:4]
fp_direct_file=dir*"15mvnm_ch1.mat"
ld_direct_file=dir*"15mvnm_ch2.mat"


resistor_value=10/5 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)

