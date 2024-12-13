include("gas_experiment_functions.jl")
direct_path="data/20241211/direct_long.mat"
direct_path="data/20241211/direct_long.txt"
L_13=load(direct_path[1:end-3] * "_L.jld2") 
plot(L_13["1"][1:100:end])