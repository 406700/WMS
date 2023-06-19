#all units in nm
using Plots
using Statistics
plotlyjs()

function plot_lor_model(amplitude)
    Γ=0.218
    x0=1547
    x=range(1546.5,step=1e-3,stop=1547.5)
    lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
    plot(x,lor_model.(x),label="model")
end
function plot_lor_prime_model(amplitude)
    Γ=0.218*1e-9
    x0=1547*1e-9
    x=range(1546.5,step=1e-3,stop=1547.5).*1e-9
    lor_prime_model(x)= amplitude*-(16*(x-x0)*Γ)/(pi*(4*(x-x0)^2+(Γ)^2)^2)/(2/(pi*Γ))*dlambda_dnu(x)
    plot(x,lor_prime_model.(x),label="model")
end

function dlambda_dnu(λ)
    nu=λ^2/3e8
end

function h(t)
    #create the chirp function
    #nb function needs to be changed with the length of t
    if t<=100
        h=log10(t)-1
    else t>100
        h=log10(t-100)-1
    end
    
end

function m(t)
    #create a square wavelength
    #nb needs to match length t
    if t<=100
        I=1
    else
        I=-1
    end
    return I
end


function simulate_fp_transmission(amplitude)
    fp_transmission_t_λ=zeros(length(t),length(center_wavelength))
    Γ=0.218
    x0=1547
    
    local lor_model(x)= amplitude*(1/pi*(1/2*Γ)/((x-x0)^2+(1/2*Γ)^2))/(2/(pi*Γ))
    for i in 1:length(center_wavelength)
        λ=λ_t.+(center_wavelength[i]) #centerered on λ0
        fp_transmission_t_λ[:,i]=I_t.*lor_model.(λ)
    end
    return fp_transmission_t_λ
end

function average_values(fp_transmission)
    FP=fp_transmission
    averages=Int(length(t)/2)
    T=averages #period over 2
    i_plus=zeros(Int(length(center_wavelength)),averages)
    i_minus=zeros(Int(length(center_wavelength)),averages)
    for i in 1:averages
        ts=(t[i])
        i_plus[:,i]=FP[ts,:].+FP[ts+T,:]
        i_minus[:,i]=FP[ts,:].-FP[ts+T,:]
    end
    i_plus=mean(i_plus,dims=2)
    i_minus=mean(i_minus,dims=2)
    return i_plus/I0,i_minus/I0
end

function L0(i_p,i_m)
     1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m)
    1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end


################################################################ 
amplitude=0.31 #determines the amplitude of the resonance, and thus it's derivative
center_wavelength=range(1546.5,step=0.01,stop=1547.5)
Δλ=0.1/2 #total chirp over 2
t=range(1,step=1,stop=200)#nb if changed, need to change functions It lambda_t
λ_t= Δλ*m.(t).*h.(t)
λ0=1547
Δν_h=2.99e8*Δλ*1e-9/(λ0*1e-9)^2*(mean(h.(t)))

ΔI=1
I0=2
I_t=ΔI*m.(t).+I0
Δim=ΔI/I0
###############################################################

fp_transmission_t_λ=simulate_fp_transmission(amplitude)
i_plus,i_minus=average_values(fp_transmission_t_λ)

L=L0.(i_plus,i_minus)
L_prime=L0_prime.(i_plus,i_minus)

#############################################################
#plot(fp_transmission_t_λ[:,1])
##
plot_lor_model(amplitude)
plot!(center_wavelength,L)
plot_lor_prime_model(amplitude)
plot!(center_wavelength*1e-9,L_prime)
dont_plot="plot = false"

