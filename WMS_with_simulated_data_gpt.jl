# Import necessary packages for plotting and statistical functions
# include("gas_experiment_functions.jl") # L L0 functions
using Interpolations,Plots,Statistics,Interpolations,MAT
# Function to plot the Lorentzian model for resonance
function plot_lor_model(amplitude)
    Γ=0.218  # Define the resonance width (FWHM)
    x0=1547  # Define the resonance center frequency
    # Generate a range of frequencies around the center frequency
    x=range(1546.5,step=1e-3,stop=1547.5)
    # Define the Lorentzian model as a function of frequency
    lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
    # Plot the Lorentzian model
    plot(x,lor_model.(x),label="model")
end

# Function to plot the derivative of the Lorentzian model
function plot_lor_prime_model(amplitude)
    Γ=0.218*1e-9  # Resonance width in meters (conversion from nm)
    x0=1547*1e-9  # Center frequency in meters
    # Generate a range of frequencies in meters
    x=range(1546.5,step=1e-3,stop=1547.5).*1e-9
    # Define the derivative of the Lorentzian model as a function of frequency
    lor_prime_model(x)= amplitude*-(16*(x-x0)*Γ)/(pi*(4*(x-x0)^2+(Γ)^2)^2)/(2/(pi*Γ))*dlambda_dnu(x)
    # Plot the derivative of the Lorentzian model
    plot(x*1e9,lor_prime_model.(x),label="model")
end

# Function to convert wavelength to frequency
function dlambda_dnu(λ)
    nu=λ^2/3e8  # Convert wavelength to frequency
end

# Function to generate the frequency chirp function h(t)
function h(t)
    # Define the chirp function with a logarithmic dependence
    if t<=100
        h=log10(t)-1
    else
        h=log10(t-100)-1
    end
end

# Function to generate a square wave function m(t) for modulation
function m(t)
    # Square wave modulation with period 200 (100 high, 100 low)
    if t<=100
        I=1
    else
        I=-1
    end
    return I
end

# Function to simulate the fiber Fabry-Perot transmission based on Lorentzian model
# function simulate_fp_transmission(amplitude)
#     fp_transmission_t_λ=zeros(length(t),length(center_wavelength))
#     Γ=0.218  # Resonance width
#     x0=1547  # Center frequency
    
#     # Define the Lorentzian model for use within this function
#     local lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
#     for i in 1:length(center_wavelength)
#         λ=λ_t.+(center_wavelength[i]) # Center on λ0 with additional modulation
#         fp_transmission_t_λ[:,i]=I_t.*lor_model.(λ)  # Calculate the transmission over time for each wavelength
#     end
#     return fp_transmission_t_λ
# end
function interpolate_hitran()
    hitran_dict=matread("hitran.mat")
    hitran_wavelength=reverse(+1 ./ hitran_dict["wavenumber"]*1e7)
    itp=LinearInterpolation(hitran_wavelength,reverse(hitran_dict["transmittance"]))
    hitran_dict=nothing
    GC.gc()
    start=minimum(hitran_wavelength)
    stop=maximum(hitran_wavelength)

    hitran_dict=matread("hitran_derivative.mat")
    itp2=LinearInterpolation(hitran_wavelength,reverse(hitran_dict["derivative_transmittance"]))

    return itp,itp2,start,stop
end
function simulate_fp_transmission(transmittance,start,stop,t,λ_t,I_t)
    data_length=1000
    wavelength=range(start=start+2*Δλ,stop=stop-2*Δλ,length=data_length)#try to avoid going of range.
    fp_transmission_t_λ=zeros(length(t),data_length)
            
        # # Define the Lorentzian model for use within this function
        # local lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
        for (index,current_wavelength) in enumerate(wavelength)
            λ=λ_t.+current_wavelength # Center on λ0 with additional modulation
            fp_transmission_t_λ[:,index]=I_t.*transmittance.(λ)  # Calculate the transmission over time for each wavelength
        end
        return fp_transmission_t_λ
end
# Function to calculate average values for sum and difference signals
function average_values(fp_transmission,data_length)
    FP=fp_transmission
    averages=Int(length(t)/2)
    T=averages # Half period of modulation
    i_plus=zeros(data_length,averages)
    i_minus=zeros(data_length,averages)
    for i in 1:averages
        ts=(t[i])
        # Calculate sum and difference signals for each wavelength
        i_plus[:,i]=FP[ts,:].+FP[ts+T,:]
        i_minus[:,i]=FP[ts,:].-FP[ts+T,:]
    end
    # Average the sum and difference signals over time
    i_plus=mean(i_plus,dims=2)
    i_minus=mean(i_minus,dims=2)
    return i_plus/I0,i_minus/I0
end

# Functions to normalize and calculate L0 and L0 prime based on average values
function L0(i_p,i_m)
     1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m)
    1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

# using Random  # Import the Random module for generating random numbers

# function simulate_fp_transmission_w_noise(amplitude, noise_level_λ, noise_level_I)
#     # The 'noise_level_λ' and 'noise_level_I' parameters determine the standard deviation of the noise
#     fp_transmission_t_λ_noise=zeros(length(t),length(center_wavelength))
#     Γ=0.218  # Resonance width
#     x0=1547  # Center frequency

#     # Define the Lorentzian model within this function
#     local lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
    
#     # Generate noisy wavelength and intensity data
#     λ_noise = Δλ .* m.(t) .* h.(t) .+ noise_level_λ .* randn(length(t))
#     I_noise = ΔI .* m.(t) .+ I0 .+ noise_level_I .* randn(length(t))
    
#     for i in 1:length(center_wavelength)
#         # Apply the noisy wavelength to the center wavelength
#         λ=λ_noise .+ center_wavelength[i]
#         # Apply the noisy intensity to the Lorentzian model
#         fp_transmission_t_λ_noise[:,i]=I_noise .* lor_model.(λ)
#     end
    
#     return fp_transmission_t_λ_noise
# end

################################################################ 
#Hitran data 
# # Main simulation parameters
# amplitude=0.57#  # Amplitude of the resonance
# center_wavelength=range(1546.5,step=0.01,stop=1547.5)
# matfile=matopen
Δλ=0.05/2  # Half of the total frequency chirp
t=range(1,step=1,stop=200)  # Time range for simulation
λ_t= Δλ*m.(t).*h.(t)  # Modulated wavelength over time
λ0=1547  # Central wavelength
Δν_h=2.99e8*Δλ*1e-9/(λ0*1e-9)^2*(mean(h.(t)))  # Normalized frequency chirp average

ΔI=1  # Amplitude modulation index
I0=4  # Average laser power
I_t=ΔI*m.(t).+I0  # Total intensity modulation over time
Δim=ΔI/I0  # Normalized amplitude modulation index
###############################################################

# Run the simulation of the FP transmission and calculate the average values
hitran,hitran_deriv,start,stop=interpolate_hitran()
fp_transmission_t_λ=simulate_fp_transmission(hitran,start,stop,t,λ_t,I_t)
i_plus,i_minus=average_values(fp_transmission_t_λ,1000) #NB data length, tmp solution

# Normalize the simulation results
L=L0.(i_plus,i_minus)
L_prime=L0_prime.(i_plus,i_minus)

plot(fp_transmission_t_λ[:,1])
wavelength=range(start=start+2*Δλ,stop=stop-2*Δλ,length=1000)#try to avoid going of range.

# Plot the Lorentzian model and the simulation results
plot(wavelength,hitran.(wavelength))
plot!(wavelength,L)
plot(wavelength,-hitran_deriv.(wavelength)/3e10)
plot!(wavelength,L_prime)

