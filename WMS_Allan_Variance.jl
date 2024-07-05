include("wms_main_functions.jl")
using Printf
rkey=Dict(zip(["λ", "ν", "L", "L_prime", "L_direct_avg", "L_prime_direct","FP_array", "LD_array","scaling_factors","Δim"],[1:10;]))
modulation_values=[5,10,15,20,25,30]

ref_results=load_results("20240110__combined_results")
ref_results2=load_results("20240110__combined_results_17")

scaling_factors=[ref_results[value][rkey["scaling_factors"]] for value in modulation_values]
scaling_factors2=[ref_results2[value][rkey["scaling_factors"]] for value in modulation_values]

ref_results=empty
ref_results2=empty
GC.gc()

#results=load_results("20240114__combined_results_stab")
results=load_results("20240126_allan_26")
results2=load_results("20240126_allan_17")
@infiltrate 

#############################################################################
# written by chatgpt verified against https://tf.nist.gov/phase/Properties/four.htm
function compute_allan_variance(V, tau)
    N = length(V) #number of sample periods times 1 millisecond.
    num_intervals = floor(Int, N / tau)
    if num_intervals<2
        error("error,see code")
      end
    interval_means = [mean(V[((i - 1) * tau + 1):(i * tau)]) for i in 1:num_intervals] #y mean V(1:tau) for each interval in V
    squared_diffs = [(interval_means[i+1] - interval_means[i])^2 for i in 1:(num_intervals - 1)] #1/2(y_n+1 - y_n)^2
    allan_variance = 0.5 * mean(squared_diffs) # <>^2
    return allan_variance
end

function get_variance(local_results,scaling_factors,do_fit=true)
    #compute allan variance. The vector step size is used as the time interval. so taus values are a range from 1 to length(vector)/2
    allan_variance=Dict()
    for (index,value) in enumerate(modulation_values)
        L_prime=local_results[value][rkey["L_prime"]]*scaling_factors[index]
        taus=1:floor(Int,length(L_prime)/2)
        if do_fit==false
            allan_variance[value]=[compute_allan_variance(L_prime,tau) for tau in taus]
        elseif do_fit==true
            fit=fitlinear(1:length(L_prime),L_prime) #
            y=L_prime.-fit.y.+fit.b
            allan_variance[value]=[compute_allan_variance(y,tau) for tau in taus]
        end
    end
    return allan_variance
end

function custom_formatter(x)
    return @sprintf("%g", x)  # This will use general number format
end

unfit_variance=get_variance(results,scaling_factors,false)
variance=get_variance(results,scaling_factors,true)

unfit_variance2=get_variance(results2,scaling_factors2,false)
variance2=get_variance(results2,scaling_factors2,true)


results=empty
results2=empty
GC.gc()

#stop=length(variance[10])
stop=250 #manually set the averaging. Confidence becomes very low. 
xticks = ([1, 10, 100], ["1", "10", "100"])
yticks=([.01,0.05],["0.01","0.05"])
pf=plot(frame=:box,title="26 mA",xticks=xticks,yticks=yticks,label="Allan Deviation", xlabel="", ylabel="",  xscale=:log10, yscale=:log10,legend=:topright,linestyle=:solid, linealpha=0.5, linewidth=4*upscale)
for (index,value) in enumerate(modulation_values)
    plot!(pf, [1:stop],sqrt.(variance[value][1:stop])*1e12,label="$value mv")  #NB sqrt 
end
display(pf)
#savefig("./figures/"*"alan_deviation_fit")

stop=250 #manually set the averaging. Confidence becomes very low. 
xticks = ([1, 10, 100], ["1", "10", "100"])
yticks=([.01,0.05],["0.01","0.05"])
p=plot(frame=:box,title="26 mA",xticks=xticks,yticks=yticks,label="Allan Deviation", xlabel="", ylabel="",  xscale=:log10, yscale=:log10,legend=:topright,linestyle=:solid, linealpha=0.5, linewidth=4*upscale)
for (index,value) in enumerate(modulation_values)

    plot!(p, [1:stop],sqrt.(unfit_variance[value][1:250])*1e12,label="$value mv")  #NB sqrt 
end
display(p)
#savefig("./figures/"*"alan_deviation_unfit")

stop=250 #manually set the averaging. Confidence becomes very low. 
xticks = ([1, 10, 100], ["1", "10", "100"])
yticks=([.01,0.05],["0.01","0.05"])
p2f=plot(frame=:box,title="17 mA",legend=false,xticks=xticks,yticks=yticks,label="Allan Deviation", xlabel="Averaging time τ (mS) ", ylabel="Allan deviation (1/THz)",  xscale=:log10, yscale=:log10,linestyle=:solid, linealpha=0.5, linewidth=4*upscale)
for (index,value) in enumerate(modulation_values)
    plot!(p2f, [1:stop],sqrt.(variance2[value][1:stop])*1e12,label="$value mv")  #NB sqrt 
end
display(p2f)
#savefig("./figures/"*"alan_deviation_fit_17")

stop=250 #manually set the averaging. Confidence becomes very low. 
xticks = ([1, 10, 100], ["1", "10", "100"])
yticks=([.01,0.05],["0.01","0.05"])
p2=plot(frame=:box,title="17 mA",legend=false,xticks=xticks,yticks=yticks,label="Allan Deviation", xlabel="Averaging time τ (mS) ", ylabel="Allan deviation (1/THz)",  xscale=:log10, yscale=:log10,linestyle=:solid, linealpha=0.5, linewidth=4*upscale)
for (index,value) in enumerate(modulation_values)
        plot!(p2, [1:250],sqrt.(unfit_variance2[value][1:250])*1e12,label="$value mv")  #NB sqrt 
end
display(p2)
#savefig("./figures/"*"alan_deviation_unfit_17")

#################################################################################################create 2x1 plot
gr()
combined_plot=plot(p,p2,layout=(2,1),link=:x,ylabel="Allan deviation (1/THz)")
combined_plot_fit=plot(pf,p2f,layout=(2,1),link=:x,ylabel="Allan deviation (1/THz)")
display(combined_plot)
savefig("./figures/"*"combined_alan_deviation.svg")

#################################################################################find the slope of the derivative function 

results=load_results("20240110__combined_results")
#results2=load_results("20240110__combined_results_17")

idx=30
L_direct=results[idx][rkey["L_direct_avg"]]
λ=results[idx][1]
scale_factor=results[idx][rkey["scaling_factors"]]

slope=calculate_derivative_slope(λ,L_direct) #dI/

λ_0=1546.6e-9
c=2.99e8
dI_dnu_d_lambda=6e-11*c/λ_0^2/140e-12
####################################################################################################calculate the resolution

deviation=[sqrt(variance[value][10]) for value in modulation_values]  ##select the averaging time here index 10 is 10 ms.
deviation2=[sqrt(variance2[value][10]) for value in modulation_values]

resolution=@. abs(deviation*c/λ_0^2/slope*1e15) #dI/dν*dν/dλ/ dI/dλ^2
resolution2=@. abs(deviation2*c/λ_0^2/slope*1e15)

deviation_unfit=[sqrt(unfit_variance[value][10]) for value in modulation_values]
deviation2_unfit=[sqrt(unfit_variance2[value][10]) for value in modulation_values]

resolution_unfit=@. abs(deviation_unfit*c/λ_0^2/slope*1e15) #dI/dν*dν/dλ/ dI/dλ^2
resolution2_unfit=@. abs(deviation2_unfit*c/λ_0^2/slope*1e15)

scaling=1/0.8
fntsm = Plots.font("sans-serif", pointsize=round(12.0*upscale)*scaling)
fntlg = Plots.font("sans-serif", pointsize=round(16.0*upscale)*scaling)
p=plot(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm*0.8,frame=:box,xlabel="modulation values (mV)",ylabel="resolution (fm)",legend=:topright,linestyle=:solid, linealpha=0.5, linewidth=4*upscale)
plot!(p,modulation_values,resolution,label="26 mA")
plot!(p,modulation_values,resolution2,label="17 mA")
savefig("./figures/"*"resolution")

p=plot(ylims=(0,70),titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, frame=:box,xlabel="modulation values (mV)",ylabel="resolution (fm)",linestyle=:solid, linealpha=0.5, linewidth=5*upscale,legend=:topright)
plot!(p,modulation_values,resolution_unfit,label="26 mA")
plot!(p,modulation_values,resolution2_unfit,label="17 mA")
savefig("./figures/"*"resolution_unfit.svg")

##############################################################################check the scaled derivatives, are they at the same transmission value?
 plotlyjs()
 p=plot()
 for (index,value) in enumerate(modulation_values)
     plot!(p,results[30][1],scaling_factors[index]*results[value][rkey["L_prime"]])
 end
 display(p)

 p=plot()
 for (index,value) in enumerate(modulation_values)
     plot!(p,results[value][1],results[value][rkey["L"]])
 end
 display(p)

# p=plot()
# for (index,value) in enumerate(modulation_values)
#     plot!(p,results[value][1],results[value][rkey["L"]])
# end
# display(p)
# # plot(results[15][1],results[15][4])
# # plot(λ,-100*λ.-λ[500])


# ######################################################################check the scaling. 
# p=plot()
# for (index,value) in enumerate(modulation_values)
#     plot!(p,results2[value][rkey["L_prime"]].*scaling_factors2[index])
# end
# display(p)
# plot!(p,plot(results[30][rkey["L_prime_direct"]]))