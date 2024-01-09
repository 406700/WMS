using DelimitedFiles, CSV, DataFrames, Plots, Statistics, MAT, LsqFit,Serialization, RollingFunctions
plotlyjs()
include("wms_main_functions.jl")


modulation_values=20#[5,10,15,20,25,30]

results=Dict()
for i in 1:length(modulation_values)
    index=modulation_values[i]
    string="data_config_"*"$index"*"ma.jl"
    include(string)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file)
    λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
    results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array) 
end

# results=Dict()
# for i in 1:length(modulation_values)
#     index=modulation_values[i]
#     string="data_config_17_"*"$index"*"mv.jl"
#     include(string)
#     #λ,ν,L,L_prime,L_direct_avg,L_prime_direct= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file)
#     λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array= load_and_process_data(fp_file,ld_file,fp_direct_file,ld_direct_file);
#     results[index]=( λ,ν,L,L_prime,L_direct_avg,L_prime_direct,FP_array, LD_array) 
# end


# Example usage
#results = load_results("results.dat")
#save_results(results,"20240106_combined_results_17ma")
#results=load_results("20240106_combined_results_17ma")

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