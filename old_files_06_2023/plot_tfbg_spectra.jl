using LsqFit
using DelimitedFiles
using Plots
plotlyjs()
cd("Data/20220722/")

function TFBG_spectrum_setup(peak_start,filename)
    #read the grating spectrum, cut it to 2 nm around the starting wavelength, normalize the percentage of loss to the starting wavelength, interpolate the results and return new vectors with 1pm resolution.
    data_OSA=readdlm(filename,',',skipstart=79)[1:end-1,:] #NB OSA has option to export with ; not ,
    grating_wavelength=data_OSA[:,1] #shorten grating spectrum to include a few peaks around λ0
    grating_spectrum=data_OSA[:,2]
    #cut spectrum to single peak. Manually adjusted.
    spectrum_start_index=findfirst(x->x>(peak_start),grating_wavelength)
    spectrum_stop_index=findfirst(x->x>(peak_start+0.5),grating_wavelength)
    minIndex=findmin(grating_spectrum[spectrum_start_index:spectrum_stop_index])[2]
    minIndex=spectrum_start_index+minIndex-1
    width=12
    grating_spectrum=grating_spectrum[minIndex-width:minIndex+width]
    max=findmax(grating_spectrum)[1]
    grating_spectrum=grating_spectrum/max
    grating_wavelength=grating_wavelength[minIndex-width:minIndex+width]
    fit=perform_fit(grating_wavelength,grating_spectrum)
    println(fit.converged)
    return grating_wavelength,grating_spectrum,fit.param
end

function perform_fit(wavelength,spectrum)
    #p=[σ,u]
    p0=[0.3,1547.4,1]

    fit = LsqFit.curve_fit(gaussian, wavelength, spectrum,p0)
    return fit
end

@. gaussian(λ, p) = 1-p[3]*1/(sqrt(2*pi)*p[1])*exp(-(1/2*(λ-p[2])^2)/p[1]^2)
files=readdir()[6:end]

gaussian_parameters=Vector{Array{Float64, 1}}(undef,length(files))

plt=plot()
for i in 1:length(files)
    wavelength_range,grating_spectrum,gaussian_parameters[i]=TFBG_spectrum_setup(1547.1,files[i])
    gui(plot!(plt,wavelength_range,grating_spectrum,label=files[i]))
    wavelength=range(wavelength_range[1],stop=wavelength_range[end],length=1000)
    plot!(wavelength,gaussian(wavelength,gaussian_parameters[i]))
end
