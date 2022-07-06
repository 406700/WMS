using LsqFit
using CurveFit
using Plots
using DelimitedFiles
using Interpolations
using EasyFit

plotlyjs()

#laser chirp properties, needst to correspond to files.
T=37.3
global const λ0=1548.252
global const ΔI=18 #ouellette
A=13.5

function perform_fit(ΔI, λ0, spectrum, volt_grating, volt_nograting, time_axis)
    function grating_extrapolate(λ)
        if λ < 1548.003047479496
            return 4
        elseif λ > 1548.59021
            return 4
        else
            return spectrum(λ)
        end
    end

    t=time_axis
    #t = range(0, 10, length = 6)#microsecond
    #p = [A, B, C, D, τ1, τ2]
    #coefficient from ouellette
    p0 = [1.4e-3, 9.79e-3, 1.02e-3, 0.0, 1.9, 15.9] #nm, mA, us
    f(t,p)=λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t)
    #
    @. model(t, p) = grating_extrapolate(λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t))/volt_nograting(t)
    fit = LsqFit.curve_fit(model, time_axis, volt_grating, p0)

    # mine(t)=f(t,fit.param)#/volt_nograting(t)
    # original(t)=f(t,p0)
    # plt=plot(time_axis,mine.(t))
    # plot!(time_axis,original.(t))

    plt2=(time,axis)
    #gui(plt)

    return fit,model
end

function get_grating_spectrum(λ0)
    data_OSA = readdlm("./Data/spec_file.csv", ',', header = true)[1]
    grating_wavelength = data_OSA[:, 1]
    spectrum=data_OSA[:, 2]
    #low_state_wavelength = findfirst(x->abs(x-λ0) <=resolution/2,grating_wavelength)
    grating_spectrum_itp=Interpolations.LinearInterpolation(grating_wavelength,spectrum)
    #grating_spectrum_itp=extrapolate(grating_spectrum_itp, 5)
    loss_at_λ0=grating_spectrum_itp(λ0)
    grating_spectrum=x->grating_spectrum_itp(x)/loss_at_λ0 #nomalize spectrum to initial wavelength

    # p0=[1.25,1548,1.0]
    # @. model(x,p) = p[1]*exp(-(x-p[2])^2/(2*p[3]^2))
    # fit = LsqFit.curve_fit(model,grating_wavelength,grating_loss,p0)
    # grating_loss_fit(t)=model(t,fit.param)
    return grating_spectrum, grating_wavelength
end

#get the oscilloscope voltage data and cut it from the trigger to the end of the pulse.
#pulse length sets the number of time steps included after the trigger
function get_oscilloscope_data()
    pulselength=25000
    volt_nograting=readdlm("./Data/373C_POWER2-2.csv",',')
    time_axis=volt_nograting[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis) #find trigger
    time_axis=time_axis[t_0:t_0+pulselength] #set pulse length

    #read the rise time of laser and detector combo. interpolate, and then curve fit, as reverse combo didn't work.
    #could use a polynomial fit to have the exact equation. Maybe look if this can be done with spline.
    volt_nograting=volt_nograting[t_0:t_0+pulselength,2]
    volt_nograting=fitspline(time_axis,volt_nograting)
    volt_nograting=LinearInterpolation(volt_nograting.x,volt_nograting.y) #x axis is shorter than original time_axis.
    #volt_nograting=LinearInterpolation(time_axis,volt_nograting)
    volt_grating=readdlm("./Data/373C_3kHz_TFBG2-2.csv",',')
    volt_grating=volt_grating[t_0:t_0+pulselength,3]


    return time_axis, volt_nograting, volt_grating
end


grating_spectrum, grating_wavelength=get_grating_spectrum(λ0)
time_axis, volt_nograting, volt_grating=get_oscilloscope_data()
int=findfirst(x->x>7.92,time_axis)
int=int-1
fit,model=perform_fit(ΔI, λ0, grating_spectrum, volt_grating[1:int], volt_nograting,time_axis[1:int])



##plotting results
# p=fit.param
# fitted_voltage=t->model(t,p)
# pfit(t)=(λ0 + ΔI*(p[1]+p[2]+(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t))/volt_nograting(t)
# plot(time_axis,fitted_voltage.(time_axis))
# plot(time_axis,volt_grating)
#
# scatter(time_axis)
#
# plot(time_axis,voltage_nograting.(t))
#

#notes:the voltage response of the laser photodiode combo (volt_nograting) must be de-noised (curve_fit)
