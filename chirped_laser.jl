using Statistics
using Plots
using DelimitedFiles
using LsqFit
using EasyFit
using Interpolations

plotlyjs()
include("chirped_laser_data.jl")
#find chirp wavelength as a  function of t, by matching the loss to the known grating profile

function calibrate_wavelength(grating_wavelength,grating_spectrum,total_loss,initial_wavelength_index)
    λ_of_t=Vector{Float64}(undef,length(total_loss))
    current_index=initial_wavelength_index
    outliers=zeros(length(total_loss))

    for i = 1:length(total_loss)
         try
            current_index = find_wavelength(total_loss[i], grating_spectrum, current_index)
            λ_of_t[i] = grating_wavelength[current_index]
         catch
             λ_of_t[i]=0
             outliers[i] = total_loss[i]
         end
    end

    return λ_of_t,outliers
end

#compares the attenuation of the single passed through the TFBG to that of the known grating spectrum. The resolution should be comparable to the noise and should ensure a steady increase in wavelength
#The function avoids jumps in wavelength by searching around the last wavelength value.

function find_wavelength(total_loss, grating_spectrum, current_index)
    #resolution = 0.05 # related to the wavelength resolution variable in TFBG_spectrum_setup
    for i = 1:length(grating_spectrum[current_index-200:end])# NB subtration starts a little below the currunt wavelength value in case of noise or reverse in direction. Will break if spectrum starts near the initial wavelength.
        if abs(grating_spectrum[i] - total_loss) < resolution
            return i
        end
    end
    #return (findmin(grating_spectrum[(current_index-200):(current_index+200)]-total_loss)[2]+current_index-1)

    #error("measured loss and spectrum loss do not match. Check resolution, λ0, or consider if grating has shifted.")
end

#osc_file= photodetector output with chirp passing through grating.
#spec file, is the grating data from the OSA

function initialize(experiment)
    # even with same oscilloscope settings trigger may not occur at exactly the same place.
    # need to read time axis from both files.
    osa_file=experiment.OSA_file
    tfbg_file=experiment.tfbg_file
    system_risetime_file=experiment.system_risetime_file
    λ0=experiment.λ0

    pulselength=Int(1e6)#1.04
    data=readdlm(tfbg_file,',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    #time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_grating=data[t_0:t_0+pulselength,3] #NB double check # of channels

    data=readdlm(system_risetime_file,',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger. Trigger should be set to volt generator.
    time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_nograting=data[t_0:t_0+pulselength,3] #NB double check files

    #account for the different losses between the two spectrum. .
    v0_grating=mean(volt_grating[1:1000])
    v0_nograting=mean(volt_nograting[1:1000])
    connectorloss=v0_grating/v0_nograting
    volt_nograting=volt_nograting.*connectorloss
    #gui(plot(time_axis,hcat(volt_nograting,volt_grating)))
    total_loss=volt_grating./volt_nograting
    #gui(plot(total_loss))
    grating_wavelength,grating_spectrum,initial_wavelength_index=TFBG_spectrum_setup(λ0,osa_file) #nb fwhm
    #gui(plot(grating_wavelength,grating_spectrum))
    result,outliers=calibrate_wavelength(grating_wavelength,grating_spectrum,total_loss,initial_wavelength_index)
    #gui(plot(time_axis,hcat(outliers,total_loss),title="outliers"))
    return result, time_axis#,fwhm

end

function TFBG_spectrum_setup(λ0,filename)
    #read the grating spectrum, cut it to 2 nm around the starting wavelength, normalize the percentage of loss to the starting wavelength, interpolate the results and return new vectors with 1pm resolution.
    data_OSA=readdlm(filename,',',skipstart=79)[1:end-1,:]#,header=true)[1]
    grating_wavelength=data_OSA[:,1] #shorten grating spectrum to include a few peaks around λ0

    spectrum_width=2 #nm
    spectrum_start_index=findfirst(x->x>(λ0-spectrum_width/2),grating_wavelength)
    spectrum_stop_index=findfirst(x->x>(λ0+spectrum_width/2),grating_wavelength)

    grating_wavelength=grating_wavelength[spectrum_start_index:spectrum_stop_index]
    grating_spectrum=data_OSA[:,2][spectrum_start_index:spectrum_stop_index]
    #gui(plot(grating_wavelength,grating_spectrum))
    #interpolate the TFBG spectrum and create new vectors with 1 pm resolution. Find the interpolated value at the exact starting wavelength for normalization
    grating_spectrum_itp=Interpolations.LinearInterpolation(grating_wavelength,grating_spectrum) #interpolate spectrum. not sure linear is correct.
    normalization_constant=grating_spectrum_itp(λ0)
    #spectral_resolution=5e-5 #wavelength step size in nm
    wavelength_range=range(grating_wavelength[1],stop=grating_wavelength[end],step=spectral_resolution) #new wavelength vector
    grating_spectrum=grating_spectrum_itp.(wavelength_range)  #new spectrum vector from wavelength range
    grating_spectrum=grating_spectrum/normalization_constant

    #get the new starting index for λ0
    initial_wavelength_index_itp=findfirst(x->x>λ0,wavelength_range)
    #fwhm=find_fwhm(wavelength_range,grating_spectrum,initial_wavelength_index_itp)
    return wavelength_range,grating_spectrum,initial_wavelength_index_itp#,fwhm
end

function find_fwhm(wavelength,grating_spectrum,initial_wavelength_index)
    spectrum_width=Int(round((0.400/spectral_resolution)))# 400pm*/resolution=number of points
    Δ=grating_spectrum[initial_wavelength_index]-grating_spectrum[initial_wavelength_index+1] #representative change in power
    search_index_1=Int(initial_wavelength_index+round(spectrum_width/2))
    search_index_2=Int(initial_wavelength_index+spectrum_width)
    stop_index=Int(findfirst(x->x-1<Δ,grating_spectrum[search_index_1:search_index_2]))
    stop_index=Int(search_index_1+stop_index-1)

    grating_spectrum=grating_spectrum[initial_wavelength_index:stop_index]
    min_index=findmin(grating_spectrum)[2] #inde
    depth=grating_spectrum[1]-grating_spectrum[min_index]
    width1=findfirst(x->abs(x-depth/2)<Δ,grating_spectrum)
    width2=findlast(x->abs(x-depth/2)<Δ,grating_spectrum)
    fwhm=abs(wavelength[width1]-wavelength[width2]) #spacing is linear so not worrying about exact indexing
    println("fwhm =$fwhm in nm")

    return fwhm
end

function perform_fit(ΔI, λ0, chirp, time_axis)
    #t=time_axis.*1e3#convert to microseconds
    p0 = [1.4, 9.79, 1.02, 0.0, 1.9, 15.9] #NB pm, mA, us
    @. model(t, p) = λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t) #NB coneverts λ0 to pm
    fit = LsqFit.curve_fit(model, time_axis, chirp, p0)

    return fit
end

global const resolution=.002 # difference between grating loss and measured loss, that are considered a match. should be related to the measurement noise.
global const spectral_resolution=5e-5 #wavelength spacing of the interpolated grating profile

#p=plot()
#for i in 8:8#:length(ChirpedLaserData.experiment_list)

chirp, time_axis=initialize(ChirpedLaserData.experiment_list[1])
#plot!(p,time_axis,chirp, xlabel="time (ms)",ylabel="wavelength (nm)",series="$i", title="chirp")
#chirp_fit = fitspline(time_axis,chirp)
#plot!(p,chirp_fit.x,chirp_fit.y)

#writedlm("./Results/20220711_$names[i]",hcat(time_axis,chirp),',')


gui(plot(time_axis,(chirp.-chirp[1]),ylabel="chirp (nm)",xlabel="time uS"))
#fit=perform_fit(ΔI,0,(chirp[1:chirp_length].-1546.6),time_axis[1:chirp_length])
# plot!(time_axis[1:chirp_length],thermal_model.(time_axis[1:chirp_length]))
# plot!(time_axis[1:chirp_length],chirp_pm[1:chirp_length])
# # fit=perform_fit(ΔI,λ0,chirp,time_axis)
# p=fit.param
# plot!(time_axis,thermal_model.(time_axis))
# plot!(time_axis,chirp)

#current problems
#1: grating using linear interpolation, is not fit to gaussian

writedlm("./Results/20220713_chirp_272_18C_long",hcat(time_axis,chirp),',')
