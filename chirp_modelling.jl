using Statistics
using Plots
using DelimitedFiles
using LsqFit
using EasyFit
using Interpolations

plotlyjs()
function perform_fit(ΔI, λ0, chirp, time_axis)
    #t=time_axis.*1e3#convert to microseconds
    p0 = [1.4, 9.79, 1.02, 0.0, 1.9, 15.9] #NB pm, mA, us
    @. model(t, p) = λ0 + ΔI*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t) #NB coneverts λ0 to pm
    fit = LsqFit.curve_fit(model, time_axis, chirp, p0)

    return fit
end

function thermal_model(t)
    #us, nm,
    #p=[1.4e-3, 9.79e-3, 1.02e-3, 0.0, 1.9, 15.9]
    #λ0=1546.653
    #ouellete 9mA
    # p = [
    #     -0.052393590263504194,
    #     16.812528394501054,
    #     8.797965541376385,
    #     0.011052997309343348,
    #     3.003738776849247,
    #     22.951300494574383,
    # ]
    #ouellete 18mA
    # P = [
    #     -0.026196790200747156,
    #     8.406264266214823,
    #     4.3989828084908655,
    #     0.005526497605411446,
    #     3.003738803033219,
    #     22.951301627810935,
    # ]
    #my chirp downsampled
    # p = [
    #     5.0975683329484,
    #     13.527426000603834,
    #     4.8097504439099845,
    #     0.1369496032399537,
    #     0.44011200682932866,
    #     2.0789718872929236e-5]
    #20220713_272_18C
    # p = [
    #     4.745444391032313,
    #     13.14126460448357,
    #     4.3747920689893816,
    #     0.020075210710208447,
    #     0.5031111290357833,
    #     8.577917156094543e-6,
    # ]
    #20220713_272_18C scraped
    # p = [
    #     0.1995269205242194,
    #     5.393775219553654,
    #     1.1636603066270863,
    #     0.0021409270782827487,
    #     1.6086100795635248,
    #     25.403833177770455,
    # ]

    #20220718
    p = [
    -6.452345933684216,
    13.959920672617752,
    1.777398906584117,
    0.00359619971207761,
    1.614980443915203,
    17.399179255509587,
    ]

    0 + 18.3*(p[1]+p[2]*(1 -exp(-t / p[5])) + p[3] * (1 - exp(-t / p[6])) + p[4] * t)
end

function power_spectrum(chirp,t,λ_range)
    power=Vector{Float64}(undef,length(λ_range)-1)
    P0=1
    for i in 1:length(λ_range)-1
        λ_range[i]
        ind1=findmin(abs.(chirp.-λ_range[i]))[2] #find t value corresponding to 1st wavelength
        ind2=findmin(abs.(chirp.-λ_range[i+1]))[2]
        power[i]=t[ind2]-t[ind1]

    end
    return power
    #gui(plot(power,ylabel="power (arbitrary units)",xlabel="time (uS)"))
end
## ouellette
ouellette=readdlm("./Data/ouellette_scraped_chirp.csv",',')
plt=plot(ouellette[:,1],ouellette[:,2])
ouellette_fit=perform_fit(8.93,0,ouellette[:,2],ouellette[:,1])
p=ouellette_fit.param

plot!(ouellette[:,1],thermal_model.(ouellette[:,1]))

# ## my chirp
# myData=readdlm("./Results/20220713_chirp_272_18C_long",',')[1:310000,:]
# t=myData[:,1]
# chirp=(myData[:,2].-myData[1,2])*1000 #normalize
# t=t[1400:end].-t[1400]
# chirp=chirp[1400:end]
# plot(chirp)
# #downsample
# downsample=range(1,length(t),step=1000)
# t_downsample=t[downsample]#[t[x] for x in downsample]
# chirp_downsample=chirp[downsample]
# myFit=perform_fit(18.3,0,myData[:,2],myData[:,1])
# p=myFit.param
#
# plot(t,thermal_model.(t))
#
#
# #try from scraped data
# scrapeData=readdlm("./Results/20220713_272_18_chirp_scraped.csv",',')
# t=scrapeData[:,1]
# chirp=scrapeData[:,2]
# plot(t,chirp)
# scrapefit=perform_fit(18.3,0,chirp,t)
# p=scrapefit.param
#
# plot(t,thermal_model.(t))
#
# ##rough power plot based on the scraped data
#
# #1 algorithm to invert function.
# Δλ=5#1 pm bin
# λ_range=range(chirp[1],stop=chirp[end],step=Δλ) #new evenly spaced x (λ) axis
#
# t_range=range(t[1],stop=t[end],length=1000000)
# chirp_itp=Interpolations.LinearInterpolation(t,chirp)
# chirp_of_t=chirp_itp.(t_range)
# #plot(t_range,chirp_of_t)
#
#
#
# ##with measured chirp
# #new evenly spaced x (λ) axis
# myData=readdlm("./Results/20220713_chirp_272_18C_long",',')[1:end,:]
# plot(myData[313800:625900,2])
# t=myData[313800:625900,1]
# chirp=(myData[313800:625900,2].-myData[625900,2])*1000 #normalize
# #t=t[1400:end].-t[1400]
# #chirp=chirp[1400:end]
# Δλ=5#1 pm bin
# λ_range=range(chirp[end],stop=chirp[1],step=Δλ)
# p=power_spectrum(chirp,t,λ_range)


##20220718

myData=readdlm("./Results/20220718_chirp_exp_2022071801",',')[1:end,:]
#plot(myData[1:386000,2])
t=myData[1050:386000,1]
chirp=(myData[1050:386000,2].-myData[1050,2])*1000 #normalize
#t=t[1400:end].-t[1400]
#chirp=chirp[1400:end]
plot(t,chirp)
my_fit=perform_fit(18.3,0,chirp,t)
p=my_fit.param
plot!(t,thermal_model.(t))
