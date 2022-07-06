using Statistics
using Plots
using DelimitedFiles
using LsqFit
plotlyjs()

#find chirp wavelength as a  function of t, by matching the loss to the known grating profile

function calibrate_wavelength(grating_wavelength,grating_loss,total_loss,initial_wavelength_index)
    λ_of_t=Vector{Float64}(undef,length(total_loss))
    current_index=initial_wavelength_index
    for i in 1:length(total_loss)
        current_index=find_wavelength(total_loss[i],grating_loss,current_index)
        try
            λ_of_t[i]=grating_wavelength[current_index]
        catch
            println(i)
        end
    end
    return λ_of_t
end


#compares the attenuation of the single passed through the TFBG to that of the known grating spectrum. The resolution should be comparable to the noise and should ensure a steady increase in wavelength
#The function avoids jumps in wavelength by searching around the last wavelength value.

function find_wavelength(total_loss, grating_loss, current_index)
    resolution = 0.001
    for i = 1:length(grating_loss[current_index-2000:end])# NB subtration starts a little below the currunt wavelength value in case of noise or reverse in direction. Will break if spectrum starts near the initial wavelength.
        if abs(grating_loss[i] - total_loss) < resolution
            return i
        end
    end
    error("no matching loss found, check if the initial wavelength is correct.")
end

#osc_file= photodetector output with chirp passing through grating.
#spec file, is the grating data from the OSA


function initialize()
    # even with same oscilloscope settings trigger may not occur at exactly the same place.
    # need to read time axis from both files.
    pulselength=Int(5e5)#1.04
    data=readdlm("./Data/20220704_grating.csv",',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    #time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_grating=data[t_0:t_0+pulselength,3]

    data=readdlm("./Data/20220704_nograting.csv",',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
    volt_nograting=data[t_0:t_0+pulselength,3]

    v0_grating=mean(volt_grating[1:1000])
    v0_nograting=mean(volt_nograting[1:1000])
    connectorloss=v0_grating/v0_nograting
    volt_nograting=volt_nograting.*connectorloss
    #gui(plot(time_axis,hcat(volt_nograting,volt_grating)))
    total_loss=volt_grating./volt_nograting
    gui(plot(total_loss))

    data_OSA=readdlm("./Data/spec_file.csv",',',header=true)[1]
    grating_wavelength=data_OSA[:,1]
    initial_wavelength_index=findfirst(x->x>λ0,grating_wavelength)
    grating_loss=data_OSA[:,2]/data_OSA[:,2][initial_wavelength_index]
    gui(plot(grating_wavelength,grating_loss))
    result=calibrate_wavelength(grating_wavelength,grating_loss,total_loss,initial_wavelength_index)
    return result,time_axis
end
# function initialize()
#
#     pulselength=Int(5e5)
#     data1=readdlm("./Data/20220704_total_ALL.csv",',',skipstart=12)#time,channel1, channel,2volts seconds
#     time_axis=data1[:,1]*1e6 #convert to us
#     t_0=findfirst(x->x>0,time_axis) #find trigger
#
#     time_axis=time_axis[t_0:t_0+pulselength] #set pulse length
#     volt_grating=data1[t_0:t_0+pulselength,3]
#     volt_nograting=data1[t_0:t_0+pulselength,9]
#
#     v0_grating=mean(volt_grating[1:1000])
#     v0_nograting=mean(volt_nograting[1:1000])
#     connectorloss=v0_grating/v0_nograting
#     volt_nograting=volt_nograting.*connectorloss
#     gui(plot(time_axis,volt_nograting)
#     total_loss=volt_grating./volt_nograting
#     gui(plot!(data1[t_0:t_0+pulselength,6]))
#
#     data_OSA=readdlm("./Data/spec_file.csv",',',header=true)[1]
#     grating_wavelength=data_OSA[:,1]
#     initial_wavelength_index=findfirst(x->x>λ0,grating_wavelength)
#     grating_loss=data_OSA[:,2]/data_OSA[:,2][initial_wavelength_index]
#     gui(plot(grating_wavelength,grating_loss))
#     result=calibrate_wavelength(grating_wavelength,grating_loss,total_loss,initial_wavelength_index)
#     return result,time_axis
# end
# # function initialize()
# #     volt_grating=readdlm("./Data/373C_3kHz_TFBG2-2.csv",',')
# #
# #     voltage_of_t=readdlm("./Data/373C_3kHz_TFBG2-2.csv", ',',skipstart=1)
# #     total_loss=voltage_of_t[:,2]
# #     time_axis=voltage_of_t[:,1]
# #         #data_OSA=readdlm("./Data/spec_file.csv",',',header=true)[1]
# #     data_OSA=readdlm("./Data/spec_file.csv",',',header=true)[1]
# #     grating_wavelength=data_OSA[:,1]
# #     initial_wavelength_index=findfirst(x->x>λ0,grating_wavelength)
# #     grating_loss=data_OSA[:,2]/data_OSA[:,2][initial_wavelength_index] # grating_spectrum normalized to the starting wavelength
# #     #gui(plot(grating_wavelength,grating_loss))
# #     #gui(plot(time_axis,total_loss))
# #     result=calibrate_wavelength(grating_wavelength,grating_loss,total_loss,initial_wavelength_index)
# #      return result,time_axis
# # end

function thermal_model(t)
    #t=t*1e3
    #p=[1.4e-3, 9.79e-3, 1.02e-3, 0.0, 1.9, 15.9]
    λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t)
end

function perform_fit(ΔI, λ0, chirp, time_axis)
    #t=time_axis.*1e3#convert to microseconds

    p0 = [1.4e-3, 9.79e-3, 1.02e-3, 0.0, 1.9, 15.9] #nm, mA, us
    @. model(t, p) = λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t)
    fit = LsqFit.curve_fit(model, time_axis, chirp, p0)

    # mine(t)=f(t,fit.param)#/volt_nograting(t)
    # original(t)=f(t,p0)
    # plt=plot(time_axis,mine.(t))
    # plot!(time_axis,original.(t))

    #plt2=(time,axis)
    #gui(plt)
    return fit
end


global const λ0 = 1548.22
global const ΔI = 18.2

chirp,time_axis=initialize()
plot!(time_axis,chirp, xlabel="time (ms)",ylabel="wavelength (nm)")
# fit=perform_fit(ΔI,λ0,chirp,time_axis)
# p=fit.param
# plot!(time_axis,thermal_model.(time_axis))
# plot!(time_axis,chirp)
