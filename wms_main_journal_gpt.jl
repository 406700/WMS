include("wms_main_functions.jl")
modulation_values=[5,10,15,20,25,30]
plotlyjs()


results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="data_config_"*"$index"*"ma.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim)  
end
save_results(results,"20240110__combined_results")

results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="data_config_17_"*"$index"*"mv.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim) 
end
save_results(results,"20240110__combined_results_17")


results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="data_config_17_"*"$index"*"mv_stab.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim) 
end
save_results(results,"20240114__combined_results_17_stab")

#modulation_values=[10,20,30]
results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="data_config_"*"$index"*"_stab.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim) 
end
save_results(results,"20240114__combined_results_stab")

#modulation_values=[10,20,30]
results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="Allan_"*"$index"*"mv.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim) 
end
save_results(results,"20240126_allan_26")

results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="Allan_17_"*"$index"*"mv.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,scale_factors,Δim) 
end
save_results(results,"20240126_allan_17")

# using NamedTuples

# # Define a NamedTuple type for your results
# @namedtuple Results(λ, ν, L, L_prime, L_direct_avg, L_prime_direct, FP_array, LD_array, scale_factors, Δim)

# results_dict = Dict{Int, Results}()

# for i in 1:length(modulation_values)
#     index = modulation_values[i]
#     string = "data_config_17_$index" * "mv_stab.jl"
#     include(string)
#     λ, ν, L, L_prime, L_direct_avg, L_prime_direct, FP_array, LD_array, scale_factors, Δim = load_and_process_data(fp_file, ld_file, fp_direct_file, ld_direct_file)
    
#     results_dict[index] = Results(λ, ν, L, L_prime, L_direct_avg, L_prime_direct, FP_array, LD_array, scale_factors, Δim)
# end


#multi_plot()
# plot_max_L()


# # Assuming 'results' is already populated as per your provided code
#derivative_scaling_factors = find_scaling_factors(results)

# # Plot scaling factor against modulation index
# plot(modulation_values[[1,2,4,5]], log2.(derivative_scaling_factors[[1,2,4,5]]), title="Scaling Factor vs Modulation Index", xlabel="Modulation Index", ylabel="Scaling Factor", legend=false)

# #visual confirmation
# #derivative_scaling_visual_check(results,derivative_scaling_factors)
# # scaling manual check: Call the function to find the scaling factor 
# #scaling_factor = find_scaling_factor(results[20][4],results[20][6])

# #chirp function checks.
#  results=Dict()
#  for i in 1:length(modulation_values)
#     index=modulation_values[i]
#     string="data_config_"*"$index"*"ma.jl"
#     include(string)
#     λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array,chirp= load_and_process_data()
#     results[index]=(FP_array,L_direct_avg,λ) 
# end
# chirp_dict=Dict()
# for i in 1:1length(results)
#     index=modulation_values[i]
#     FP_array,L_direct_avg,λ=results[index]
#    chirp_dict[index]=simple_chirp_estimate(FP_array,L_direct_avg, λ,0.2)
# end

# index=modulation_values[4]
# FP_array,L_direct_avg,λ=results[index]
# @time chirp=simple_chirp_estimate(FP_array,L_direct_avg, λ,0.2)

# chirp=simple_chirp_estimate(FP_array,L_direct_avg, λ,0.2)
# plot(chirp)


# foo=exponential_fit(rollmean(FP_array[31250:end,350],100))
# plot!(rollmean(FP_array[31250:end,350],100))
# plot(chirp)

#scratch
#plot(modulation_values,derivative_scaling_factors)
#savefig("scaling_factor_v_modulation")

#multi_plot(results,scale_factor,modulation,save)
λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array=results[30]
plot(rollmean(LD_array[:,500],200)) 




#temp
# function f(τ)
#     T=1
#     1-2/(1-exp(-T/(2*τ)))+(4*exp(T/(2*τ))τ)/(T)
# end

# function f2(τ)
#     T=1
#     1-4*τ/T-2/(1-exp(-T/(2*τ)))
# end
# τ=0.5:0.01:1.5
# using Plots
# plot(τ,f2.(0.1:0.01:1))
