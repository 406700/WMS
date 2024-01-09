
dir="/home/michael/Documents/Experimental_Data/WMS_final_paper/20240104/"
fp_file =  dir*"10mv_ch1.mat"
ld_file = dir* "10mv_ch2.mat"#"path/to/first/ld/file.mat"

#FP 
FP_frequency=0.5 #hertz
M=999
#measured for .2 V modulation.
FP_min_wavelength=1546.286e-9 #NB
FP_max_wavelength=1547.078e-9 #NB
FP_shift_0= (FP_max_wavelength-FP_min_wavelength)
FP_modulation_amplitude=0.3 #volts
#adjust for modulation FP_modulation_amplitude
FP_shift=(FP_modulation_amplitude/0.2)*FP_shift_0
delta_shift=FP_shift-FP_shift_0
FP_min_wavelength=FP_min_wavelength-abs(delta_shift/2)
FP_max_wavelength=FP_max_wavelength+abs(delta_shift/2)

#LD 
c=2.99e8
average_ν=c/1546.920e-9
ν_0=c/1546.896e-9


#direct files
#f_direct=readdir(dir)[3:4]
fp_direct_file=dir*"10mvnm_ch1.mat"
ld_direct_file=dir*"10mvnm_ch2.mat"


resistor_value=10/5 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)

