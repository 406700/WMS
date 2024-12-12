# include("gas_experiment_functions.jl")
# using LsqFit


# function fit_to_hitran(x,L,itp,p0)
#     # Define the fitting model with parameters [c, chi]
#     @. model(p,x) = p[1] * (1 - p[2] * (1 - itp(x)))
#     fit = curve_fit(model,x,L,p0)

#     # Extract fitted parameters
#     fitted_params = fit.param
#     return fitted_params
# end
# ############################################################################################

# det_data_path="data/20241202/01_bar_long.mat" 
# L = load(det_data_path[1:end-3] * "_L_prime.jld2")["1"] 

# interpolated_hitran,_,start_wavelength,stop_wavelength=interpolate_hitran()
# #use manual selection of wavelength start,stop and center position 
# wavelength_scale=collect(range(start=start_wavelength,stop=1547.5,length=length(L)))#note, the hitran scan and L scan can be a different length, so choose the wavelength scaling factor based on this.
# C=1.0
# chi=0.5
# p0=[C,chi]
# fit_params=fit_to_hitran(wavelength_scale,L,interpolated_hitran,p0)

# interpolated_hitran.(wavelength_scale)
include("gas_experiment_functions.jl")
using Optim

function fit_to_hitran(x, L, itp, p0)
    # Define the model with parameters [c, chi]
    model(p, x) = p[1] * (1 .- p[2] .* (1 .- itp(x)))

    # Objective function: sum of squared residuals
    function obj(p)
        residuals = L .- model(p, x)
        return sum(residuals.^2)
    end

    # Perform optimization using, for example, LBFGS (you could try BFGS as well)
    result = optimize(obj, p0, LBFGS())

    # Extract optimized parameters
    fitted_params = Optim.minimizer(result)
    return fitted_params
end

############################################################################################

det_data_path = "data/20241202/01_bar_long.mat" 
L = load(det_data_path[1:end-3] * "_L_prime.jld2")["1"]

interpolated_hitran, _, start_wavelength, stop_wavelength = interpolate_hitran()

# Use manual selection of wavelength start, stop, and center position
wavelength_scale = collect(range(start = start_wavelength, stop = 1547.5, length = length(L)))

C = 1.0
chi = 0.5
p0 = [C, chi]

fit_params = fit_to_hitran(wavelength_scale, L, interpolated_hitran, p0)

interpolated_values = interpolated_hitran.(wavelength_scale)
fit_params, interpolated_values
