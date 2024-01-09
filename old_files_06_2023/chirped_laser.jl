using Statistics
using Plots
using DelimitedFiles
using LsqFit
using EasyFit
using Interpolations

plotlyjs()
include("chirped_laser_data.jl")
#find chirp wavelength as a  function of t, by matching the loss to the known grating profile

function calibrate_wavelength(grating_wavelength,grating_spectrum,normalized_loss,minIndex)
    λ_of_t=Vector{Float64}(undef,length(normalized_loss))
    current_index=1#initial_wavelength_index
    outliers=zeros(length(normalized_loss))
    before_minimum=true

    wait=0
    for i = 1:length(normalized_loss)
         # try
        current_index, before_minimum,wait= find_wavelength(normalized_loss[i], grating_spectrum, current_index,before_minimum,minIndex,wait)
        λ_of_t[i] = grating_wavelength[current_index]
         # catch
         #     λ_of_t[i]=0
         #     outliers[i] = normalized_loss[i]
         #     error("catch")
         # end
    end

    return λ_of_t,outliers
end

#compares the attenuation of the single passed through the TFBG to that of the known grating spectrum. The resolution should be comparable to the noise and should ensure a steady increase in wavelength
#The function avoids jumps in wavelength by searching around the last wavelength value.

function find_wavelength(normalized_loss, grating_spectrum, current_index,before_minimum,minIndex,wait)
    #resolution = 0.05 # related to the wavelength resolution variable in TFBG_spectrum_setup
    # for i = 1:length(grating_spectrum[current_index-200:end])# NB subtration starts a little below the currunt wavelength value in case of noise or reverse in direction. Will break if spectrum starts near the initial wavelength.
    #     if abs(grating_spectrum[i] - normalized_loss) < resolution
    #         return i
    #     end
    # end

    # if (abs(1-grating_spectrum[current_index])<=.016213) & (wait>2500) #wait counter ensures the minimum doesn't toggle rapidly back and fourth. If actual data switches near minimum this will cause a problem NB what value for < comparison=
    #     before_minimum=!before_minimum
    #     wait=0
    # end
    # wait+=1

    grating_spectrum=abs.(grating_spectrum.-normalized_loss)
    if before_minimum
        current_index= Int(findmin(grating_spectrum[1:minIndex])[2])
    else
        current_index= Int(findmin(grating_spectrum[minIndex:end])[2]+minIndex-1)
    end

    #NB conditional statements need care

    return current_index,before_minimum,wait
    #error("measured loss and spectrum loss do not match. Check resolution, λ0, or consider if grating has shifted.")
end

#osc_file= photodetector output with chirp passing through grating.
#spec file, is the grating data from the OSA

function initialize(powerMinimum,experiment)
    # even with same oscilloscope settings trigger may not occur at exactly the same place.
    # need to read time axis from both files.
    osa_file=experiment.OSA_file
    tfbg_file=experiment.tfbg_file
    system_risetime_file=experiment.system_risetime_file
    λ0=experiment.λ0

    pulselength=Int(4e5)#limit the time/ length of data
    data=readdlm(tfbg_file,',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    volt_grating=data[t_0:t_0+pulselength,3] #NB double check # of channels

    data=readdlm(system_risetime_file,',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger. Trigger should be set to volt generator.
    time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_nograting=data[t_0:t_0+pulselength,3] #NB double check files
    ##gui(plot(volt_nograting))
    #account for the different losses between the two spectrum. .
    #v0_grating=mean(volt_grating[1:1000])
    v0_high=mean(volt_nograting[100000:200000]) #NB .0005 V difference in value due to slope at top of pulse.  #average _  points at the top of the pulse.

    #gui(plot(hcat(volt_nograting,volt_grating)))
    normalized_loss=volt_grating./(powerMinimum/v0_high*volt_nograting)
    gui(plot(normalized_loss))
    grating_wavelength,grating_spectrum,minIndex=TFBG_spectrum_setup(λ0,osa_file) #nb fwhm
    ##gui(plot(grating_wavelength,grating_spectrum))
    result,outliers=
    ##gui(plot(time_axis,hcat(outliers,normalized_loss),title="outliers"))
    return grating_wavelength,grating_spectrum,normalized_loss,minIndex,time_axis

end
function TFBG_spectrum_setup(λ0,filename)
    #read the grating spectrum, cut it to 2 nm around the starting wavelength, normalize the percentage of loss to the starting wavelength, interpolate the results and return new vectors with 1pm resolution.
    data_OSA=readdlm(filename,',',skipstart=79)[1:end-1,:] #NB OSA has option to export with ; not ,
    grating_wavelength=data_OSA[:,1] #shorten grating spectrum to include a few peaks around λ0
    #gui(plot(grating_wavelength))
    spectrum_width=0.8 #nm
    spectrum_start_index=findfirst(x->x>(λ0-spectrum_width/4),grating_wavelength)
    spectrum_stop_index=findfirst(x->x>(λ0+spectrum_width*3/4),grating_wavelength)

    grating_wavelength=grating_wavelength[spectrum_start_index:spectrum_stop_index]
    grating_spectrum=data_OSA[:,2][spectrum_start_index:spectrum_stop_index]
    ##gui(plot(grating_wavelength,grating_spectrum))
    #interpolate the TFBG spectrum and create new vectors with 1 pm resolution. Find the interpolated value at the exact starting wavelength for normalization
    grating_spectrum_itp=Interpolations.LinearInterpolation(grating_wavelength,grating_spectrum) #interpolate spectrum. not sure linear is correct.
    #spectral_resolution=5e-5 #wavelength step size in nm
    wavelength_range=range(grating_wavelength[1],stop=grating_wavelength[end],step=spectral_resolution) #new wavelength vector
    grating_spectrum=grating_spectrum_itp.(wavelength_range)  #new spectrum vector from wavelength range
    #gui(plot(wavelength_range,grating_spectrum))
    min=findmin(grating_spectrum)
    minIndex=min[2] #index of minimum element
    grating_spectrum=grating_spectrum/min[1] #divide by value of grating minimum
    #fwhm=find_fwhm(wavelength_range,grating_spectrum,initial_wavelength_index_itp)
    #gui(plot(wavelength_range,grating_spectrum))

    return wavelength_range,grating_spectrum,minIndex#,fwhm
end

global const resolution=.002 # difference between grating loss and measured loss, that are considered a match. should be related to the measurement noise.
global const spectral_resolution=5e-5 #wavelength spacing of the interpolated grating profile

#p=plot()
#for i in 8:8#:length(ChirpedLaserData.experiment_list)
grating_wavelength,grating_spectrum,normalized_loss,minIndex,time_axis=initialize(0.03725,ChirpedLaserData.experiment_list[3])
chirp,outliers=calibrate_wavelength(grating_wavelength,grating_spectrum,normalized_loss,minIndex)
plot(chirp)
plot(normalized_loss)
#plot!(p,time_axis,chirp, xlabel="time (ms)",ylabel="wavelength (nm)",series="$i", title="chirp")
#chirp_fit = fitspline(time_axis,chirp)
#plot!(p,chirp_fit.x,chirp_fit.y)

#gui(plot(time_axis,(chirp.-chirp[1]),ylabel="chirp (nm)",xlabel="time uS"))
#fit=perform_fit(ΔI,0,(chirp[1:chirp_length].-1546.6),time_axis[1:chirp_length])
# plot!(time_axis[1:chirp_length],thermal_model.(time_axis[1:chirp_length]))
# plot!(time_axis[1:chirp_length],chirp_pm[1:chirp_length])
# # fit=perform_fit(ΔI,λ0,chirp,time_axis)
# p=fit.param
# plot!(time_axis,thermal_model.(time_axis))
# plot!(time_axis,chirp)

#current problems
#1: grating using linear interpolation, is not fit to gaussian

#writedlm("./Results/20220718_chirp_exp_2022071801",hcat(time_axis,chirp),',')
