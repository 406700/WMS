include("gas_experiment_functions.jl") #just to load packages. should edit this to only load whats needed. 
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

function get_variance(local_results,scaling_factors,do_fit=false)
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


det_data_path="data/20241202/03bar_stab.mat"

L_dict=load(det_data_path[1:end-3]*"_L.jld2")
L_prime_dict=load(det_data_path[1:end-3]*"_L_prime.jld2")


