using DelimitedFiles
using Plots
using Statistics
#cd("/home/m/OneDrive/Experimental_Data/20230620_stability")
cd("/home/m/OneDrive/Experimental_Data/WMS_final_paper/20261226_stability")
f=readdir()
L_prime_final=readdlm("20mvstab2_L_prime")
t=1:1:995

#############################################################################chatgpt code:
using Plots

function compute_allan_variance(V, tau)
    N = length(V) #number of sample periods times 1 millisecond.
    num_intervals = floor(Int, N / tau)
    if num_intervals < 2
        return NaN # Return NaN for insufficient data
    end
    interval_means = [mean(V[((i - 1) * tau + 1):(i * tau)]) for i in 1:num_intervals]
    squared_diffs = [(interval_means[i+1] - interval_means[i])^2 for i in 1:(num_intervals - 1)]
    allan_variance = 0.5 * mean(squared_diffs)
    return allan_variance
end

function plot_allan_variance(V, taus) #taus is a range of times to compute the variance for
    variances = [compute_allan_variance(V, tau) for tau in taus]

    plot(taus, variances, label="Allan Variance", xlabel="τ (Tau)", ylabel="Allan Variance", title="Allan Variance vs Tau", xscale=:log10, yscale=:log10)
end

plot_allan_variance(L_prime_final,t)
plot(L_prime_final)
floor(Int, length(L_prime_final) / 10)