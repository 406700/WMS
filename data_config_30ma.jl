
dir="/home/michael/Documents/Experimental_Data/WMS_final_paper/20240104/"
fp_file =  dir*"30mv_ch1.mat"
ld_file = dir* "30mv_ch2.mat"#"path/to/first/ld/file.mat"

#FP 
FP_frequency=0.5 #hertz
FP_shift=3.63*1e-9 #nm*volt *1e-9m/nm
FP_modulation_amplitude=0.3 #volts
FP_shift=(FP_modulation_amplitude)*FP_shift

#LD 
c=2.99e8
ν_0=c/1546.637e-9

#average_ν=c/1546.920e-9


#direct files
#f_direct=readdir(dir)[3:4]
fp_direct_file=dir*"30mvnm_ch1.mat"
ld_direct_file=dir*"30mvnm_ch2.mat" #NB data file corrupted


resistor_value=10/5 #Rch1/Rch2 Set the resistors at the channel input, which scale the voltages.
sample_interval=1.6e-8 #time[2]-time[1] #seconds
data_points_per_period=Int(0.001/sample_interval)


#measured for .2 V modulation.
# FP_min_wavelength=1546.286e-9 #NB
# FP_max_wavelength=1547.078e-9 #NB
# FP_shift_0= (FP_max_wavelength-FP_min_wavelength)
# delta_shift=FP_shift-FP_shift_0
# FP_min_wavelength=FP_min_wavelength-abs(delta_shift/2)
# FP_max_wavelength=FP_max_wavelength+abs(delta_shift/2)
#adjust for modulation FP_modulation_amplitude
