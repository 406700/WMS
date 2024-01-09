# Import necessary packages for plotting and statistical functions
using Plots
using Statistics
plotlyjs()  # Choose PlotlyJS as the backend for interactive plotting

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
function simulate_fp_transmission(amplitude)
    fp_transmission_t_λ=zeros(length(t),length(center_wavelength))
    Γ=0.218  # Resonance width
    x0=1547  # Center frequency
    
    # Define the Lorentzian model for use within this function
    local lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
    for i in 1:length(center_wavelength)
        λ=λ_t.+(center_wavelength[i]) # Center on λ0 with additional modulation
        fp_transmission_t_λ[:,i]=I_t.*lor_model.(λ)  # Calculate the transmission over time for each wavelength
    end
    return fp_transmission_t_λ
end

# Function to calculate average values for sum and difference signals
function average_values(fp_transmission)
    FP=fp_transmission
    averages=Int(length(t)/2)
    T=averages # Half period of modulation
    i_plus=zeros(Int(length(center_wavelength)),averages)
    i_minus=zeros(Int(length(center_wavelength)),averages)
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

using Random  # Import the Random module for generating random numbers

function simulate_fp_transmission_w_noise(amplitude, noise_level_λ, noise_level_I)
    # The 'noise_level_λ' and 'noise_level_I' parameters determine the standard deviation of the noise
    fp_transmission_t_λ_noise=zeros(length(t),length(center_wavelength))
    Γ=0.218  # Resonance width
    x0=1547  # Center frequency

    # Define the Lorentzian model within this function
    local lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
    
    # Generate noisy wavelength and intensity data
    λ_noise = Δλ .* m.(t) .* h.(t) .+ noise_level_λ .* randn(length(t))
    I_noise = ΔI .* m.(t) .+ I0 .+ noise_level_I .* randn(length(t))
    
    for i in 1:length(center_wavelength)
        # Apply the noisy wavelength to the center wavelength
        λ=λ_noise .+ center_wavelength[i]
        # Apply the noisy intensity to the Lorentzian model
        fp_transmission_t_λ_noise[:,i]=I_noise .* lor_model.(λ)
    end
    
    return fp_transmission_t_λ_noise
end

################################################################ 
# Main simulation parameters
amplitude=0.57#  # Amplitude of the resonance
center_wavelength=range(1546.5,step=0.01,stop=1547.5)
Δλ=0.1/2  # Half of the total frequency chirp
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
fp_transmission_t_λ=simulate_fp_transmission(amplitude)
i_plus,i_minus=average_values(fp_transmission_t_λ)

# Normalize the simulation results
L=L0.(i_plus,i_minus)
L_prime=L0_prime.(i_plus,i_minus)


############################################################
#with random noise 

# Define noise levels for wavelength and intensity
# noise_level_λ = 0.1  # Adjust the noise level as needed for wavelength
# noise_level_I = 0.5   # Adjust the noise level as needed for intensity

# # Run the noisy simulation
# fp_transmission_t_λ_noise = simulate_fp_transmission_w_noise(amplitude, noise_level_λ, noise_level_I)
# i_plus,i_minus=average_values(fp_transmission_t_λ_noise)

# # Normalize the simulation results
# L=L0.(i_plus,i_minus)
# L_prime=L0_prime.(i_plus,i_minus)
#############################################################
# Uncomment below to plot the transmission for the first center wavelength
plot(fp_transmission_t_λ[:,1])

# Plot the Lorentzian model and the simulation results
plot_lor_model(amplitude)
plot!(center_wavelength,L)
title!("Simulated")
plot_lor_prime_model(amplitude)
plot!(center_wavelength,L_prime)
title!("simulated")