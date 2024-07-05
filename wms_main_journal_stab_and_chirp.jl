include("wms_main_functions.jl")
rkey=Dict(zip(["λ", "ν", "L", "L_prime", "L_direct_avg", "L_prime_direct","FP_array", "LD_array","scale_factor","Δim"],[1:10;]))

results=load_results("20240114__combined_results_stab")
results2=load_results("20240114__combined_results_17_stab")

modulation_values=[5,10,15,20,25,30]

############################################calculate I0_direct
I_0=mean([mean(results[15][rkey["LD_array"]][:,value]) for value in 1:995])
I_02=mean([mean(results2[15][rkey["LD_array"]][:,value]) for value in 1:995])

##############################################calculate chirp
chirp_result=Dict()
for (index,value) in enumerate(modulation_values)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factor,Δim=results[value] 
    FP_average=zeros(length(FP_array[:,1]))
    for i in 1:length(FP_array[:,1])
        FP_average[i]=mean(FP_array[i,:])
    end
    chirpup,chirpdown=direct_chirp_calculation(λ,L_direct_avg,FP_average[1:31000])
    chirp_result[value]= chirpup.-chirpup[1] ,chirpdown.-chirpdown[1]  
end

chirp_result2=Dict()
for (index,value) in enumerate(modulation_values)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factor,Δim=results2[value] 
    FP_average=zeros(length(FP_array[:,1]))
    for i in 1:length(FP_array[:,1])
        FP_average[i]=mean(FP_array[i,:])
    end
    chirpup,chirpdown=direct_chirp_calculation(λ,L_direct_avg,FP_average[1:31000])
    chirp_result2[value]= chirpup.-chirpup[1] ,chirpdown.-chirpdown[1]  
end

######plot all the chirps for 26mA
time=range(start=0.0,step=0.000000016,length=length(chirp_result[5][1]))*1000
p=plot(xlabel=legend=:bottomright,frame=:box,linestyle=:solid, linealpha=0.5, linewidth=4*upscale,grid=false)
for (index,value) in enumerate(modulation_values)
    plot!(p,time,chirp_result[value][1]./maximum(abs.(chirp_result[value][1])),label="$value mV")
end
display(p)
#savefig("./figures/"*"normalized_chirp")


p=plot(ylabel="chirp (pm)",xlabel="",legend=:bottomright,frame=:box,linestyle=:solid, linealpha=0.5, linewidth=4*upscale,grid=false)
for (index,value) in enumerate(modulation_values)
    plot!(p,time,chirp_result[value][1]*1e12,label="$value mV")
end
display(p)
savefig("./figures/"*"chirp_unnormalized")

#######compare chirps betweem 17 and 26mA

p2=plot(ylabel="normalized chirp",xlabel="time (mS)",legend=:bottomright,frame=:box,linestyle=:solid, linealpha=0.5, linewidth=4*upscale)
for (index,value) in enumerate([10,20])
    plot!(p2,time,chirp_result[value][1]./maximum(abs.(chirp_result[value][1])),label="26 mA $value mV")
    plot!(p2,time,chirp_result2[value][1]./maximum(abs.(chirp_result2[value][1])),label="17mA $value mV")
end
display(p2)
savefig("./figures/"*"chirp_16_26")

combined_plot = plot(p, p2, layout=(2, 1),link=:x)
display(combined_plot)
savefig("./figures/combined_chirp_plot")

#######################create a 2x1 plot of the max chirp and Delta I

###twin plot of chirp and delta im
gr()  # Set the backend to GR
# Create the first plot with the actual data
scaling=1/0.8
fntsm = Plots.font("sans-serif", pointsize=round(12.0*upscale)*scaling)
fntlg = Plots.font("sans-serif", pointsize=round(16.0*upscale)*scaling)

p=plot(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm,xlabel="modulation voltage (mV)", ylabel="chirp (pm)",xlims=(0, maximum(modulation_values)),ylims=(0, 90),frame=:box,linestyle=:solid, linealpha=0.5, linewidth=4*upscale,grid=false, left_margin=5mm,right_margin=5mm,buttom_margin=5mm)

plot!(p,modulation_values, [chirp_result[value][1][end] for value in modulation_values]*1e12, label="chirp 26mA", color=:blue)
# Add the second data series to the same plot
plot!(p,modulation_values, [chirp_result2[value][1][end] for value in modulation_values]*1e12, label="chirp 17mA", color=:green)

# Create a secondary axis for delta but do not create a label yet
p2 = twinx()
plot!(p2, modulation_values, [results[value][end][1]*I_0 for value in modulation_values], color=:blue, linestyle=:dash, label="", ylabel="ΔI (V)",xlims=(0, maximum(modulation_values)),ylims=(0,3),frame=:box,grid=true,titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
plot!(p2, modulation_values, [results2[value][end][1]I_02 for value in modulation_values], color=:green,linestyle=:dash, label="")
# Add dummy series for the secondary axis data (for legend purposes)
plot!(p, [], [], label="ΔI 26mA", color=:blue,linestlye=:dash)
plot!(p, [], [], label="ΔI 17mA", color=:green,linestyle=:dash)


# Display the combined plot
plot(p, p2, layout=(1,1)) 
savefig("./figures/"*"chirp_vs_modulation_16_26.svg")

p=plot(ylabel="chirp (pm)",xlabel="modulation voltage",legend=:bottomright,xlims=(0, maximum(modulation_values)))
plot!(p,modulation_values,[chirp_result[value][1][end] for value in modulation_values]*1e12,label="26mA")
plot!(p,modulation_values,[chirp_result2[value][1][end] for value in modulation_values]*1e12,label="17mA")
############################################################################################################# fit the chirp 
chirp_fitted=Dict()
for (index,value) in enumerate(modulation_values)
    chirp =chirp_result[value][1]
    x=1:length(chirp)
    fit=fitlinear(x[15000:end],chirp[15000:end])
    chirp_fitted[value]=fit.y[end]
end
plot(modulation_values,[chirp_fitted[value] for value in modulation_values])

chirp_fitted_17=Dict()
for (index,value) in enumerate(modulation_values)
    chirp =chirp_result2[value][1]
    x=1:length(chirp)
    fit=fitlinear(x[15000:end],chirp[15000:end])
    chirp_fitted_17[value]=fit.y[end]
end
plot!(modulation_values,[chirp_fitted_17[value] for value in modulation_values])

############################################################################################################# create figure showing frequency noise
FP_array=results[20][rkey["FP_array"]]
FP_rollmean= mapslices(col -> rollmean(col, 200), FP_array, dims=1)
time=range(start=0.0,step=0.000000016,length=length(FP_rollmean[:,1]))*1000
p=plot( xlabel="time (mS)", ylabel=L"I_{fp}/I_{ld}",frame=:box,linestyle=:solid, linealpha=0.5, linewidth=4*upscale,grid=false)
for index in 1:5
    plot!(p,time,FP_rollmean[:,index], label="cycle $index")
end
display(p)
savefig("./figures/frequency_noise.svg")



########################################################################################calculate chirp
#chirp_dict=multi_chirp_estimate(300,500,results,modulation_values)
#λ,_,_,_,L_direct_avg,_,FP_array, LD_array=results[20] 
  
#foo=simple_chirp_estimate_stab(FP_array[:,570:581],L_direct_avg,λ)
#plot(FP_array[:,481])
#plot(LD_array[:,481] )

#scratch
# p=plot()
# for (index,value) in enumerate(modulation_values)
#     plot!(p,scaling_factors[index]*results[value][rkey["L_prime"]])
# end
# display(p)


# ############################<>#######################################################################old

    # FP_rollmean= mapslices(col -> rollmean(col, 200), FP_array, dims=1)

    # #####calculate a chirp
    # @time x=simple_chirp_estimate_stab(FP_rollmean[:,200:500],L_direct_avg,λ)
    # plot(x)
    # mean(x)
    # plot(rollmean_matrix[:,1])
    # #chirp_dict=Dict(30=>130,25=>116,20=>64,15=>63,10=>49,5=>22)
    # #recalculated after lambda redefinition (error in the loaded files definitaion of lambda) from stability measurements: 
    # plot([10,20,30],[34,64,97],legend=false)
    # plot([5,10,15,20,25,30],[22,46,63,64,116,130])

    # ##create plots for estimating the chirp
    # x=zeros(length(FP_array[:,1]))
    # for i in 1:length(FP_array[:,1])
    #     x[i]=mean(FP_array[i,:])
    # end
    # plot(x)
    # plot(FP_rollmean[:,1:5])
    # plot(rollmean(LD_array[:,1],200))

    # #pigtail measured in the lab around 75-78 cm. possible error in lambda shift. Getting 235 vs 214 pm fwhm, but not sure how precise the measurement on the OSA is.
    # #10,v 11,.42,.46,true =75.8
    # #20mv stab 15, 0.4, 0.39 false=73.9
    # #30mv stab 21 .37 .305, false  = 76
    # get_resonator_length(λ,L_direct_avg,11,.42,.46,true)
    

    # simple_chirp_estimate_stab(FP_rollmean[:,1:25],L_direct_avg,λ)

    ####################################################################################direct chirp estimate without fitting.
 =#
