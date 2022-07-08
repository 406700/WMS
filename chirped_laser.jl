using Statistics
using Plots
using DelimitedFiles
using LsqFit
using Interpolations
plotlyjs()
include("chirped_laser_data.jl")
#find chirp wavelength as a  function of t, by matching the loss to the known grating profile

function calibrate_wavelength(grating_wavelength,grating_spectrum,total_loss,initial_wavelength_index)
    λ_of_t=Vector{Float64}(undef,length(total_loss))
    current_index=initial_wavelength_index
    for i in 1:length(total_loss)
        current_index=find_wavelength(total_loss[i],grating_spectrum,current_index)
        λ_of_t[i]=grating_wavelength[current_index]
    end
    return λ_of_t
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
    error("measured loss and spectrum loss do not match. Check resolution, λ0, or consider if grating has shifted.")
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

    pulselength=Int(5e5)#1.04
    data=readdlm(tfbg_file,',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    #time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_grating=data[t_0:t_0+pulselength,3] #index three for ALL.csv files

    data=readdlm(system_risetime_file,',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_nograting=data[t_0:t_0+pulselength,2] #index 2 due to way files are saved

    #account for the different losses between the two spectrum. methodology, not code, may need review.
    v0_grating=mean(volt_grating[1:1000])
    v0_nograting=mean(volt_nograting[1:1000])
    connectorloss=v0_grating/v0_nograting
    volt_nograting=volt_nograting.*connectorloss
    gui(plot(time_axis,hcat(volt_nograting,volt_grating)))
    total_loss=volt_grating./volt_nograting
    gui(plot(total_loss))
    grating_wavelength,grating_spectrum,initial_wavelength_index=TFBG_spectrum_setup(λ0,osa_file)
    gui(plot(grating_wavelength,grating_spectrum))
    result=calibrate_wavelength(grating_wavelength,grating_spectrum,total_loss,initial_wavelength_index)

    return result,time_axis
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
    return wavelength_range,grating_spectrum,initial_wavelength_index_itp
end

function thermal_model(t)
    #us, nm,
    #p=[1.4e-3, 9.79e-3, 1.02e-3, 0.0, 1.9, 15.9]
    λ0=1546.653
    p=[]
    λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t)
end

function perform_fit(ΔI, λ0, chirp, time_axis)
    #t=time_axis.*1e3#convert to microseconds

    p0 = [1.4, 9.79, 1.02, 0.0, 1.9, 15.9] #nm, mA, us
    @. model(t, p) = λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t) #NB coneverts λ0 to pm
    fit = LsqFit.curve_fit(model, time_axis, chirp, p0)

    # mine(t)=f(t,fit.param)#/volt_nograting(t)
    # original(t)=f(t,p0)
    # plt=plot(time_axis,mine.(t))
    # plot!(time_axis,original.(t))

    #plt2=(time,axis)
    #gui(plt)
    return fit
end

global const resolution=.05 # difference between grating loss and measured loss, that are considered a match. should be related to the measurement noise.
global const spectral_resolution=5e-3 #wavelength spacing of the interpolated grating profile

p=plot()
for i in 3:3#3:length(ChirpedLaserData.experiment_list)
    experiment=ChirpedLaserData.experiment_list[i]
    chirp,time_axis=initialize(experiment)
    plot!(p,time_axis,chirp, xlabel="time (ms)",ylabel="wavelength (nm)",series="$i")
end
plot(p)


#fit=perform_fit(ΔI,0,(chirp[1:chirp_length].-1546.6),time_axis[1:chirp_length])
# plot!(time_axis[1:chirp_length],thermal_model.(time_axis[1:chirp_length]))
# plot!(time_axis[1:chirp_length],chirp_pm[1:chirp_length])
# # fit=perform_fit(ΔI,λ0,chirp,time_axis)
# p=fit.param
# plot!(time_axis,thermal_model.(time_axis))
# plot!(time_axis,chirp)

#writedlm("./Results/chirp_27_16_281.csv",hcat(time_axis,chirp),',')
