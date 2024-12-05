include("gas_experiment_functions.jl")
 #just to load packages. should edit this to only load whats needed. 
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
function compute_allan_variance_fast(V, tau)
    N = length(V)
    num_intervals = div(N, tau,RoundDown)  # Use integer division directly

    if num_intervals < 2
        error("Insufficient data: tau is too large, or the dataset is too small to compute Allan variance.")
    end

    # Precompute segment boundaries
    indices = collect(1:tau:(num_intervals * tau + 1))
    
    # Calculate interval means efficiently
    interval_means = [mean(@view V[indices[i]:(indices[i + 1] - 1)]) for i in 1:(num_intervals)]
    
    # Compute squared differences and mean
    squared_diffs = @. (interval_means[2:end] - interval_means[1:end-1])^2
    allan_variance = 0.5 * mean(squared_diffs)
    
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

L=load(det_data_path[1:end-3]*"_L.jld2")["1"]
L_prime=load(det_data_path[1:end-3]*"_L_prime.jld2")["1"]

data= @view L_prime[]#[20000:40000]
tau_s=1:floor(Int,length(data)/100)

allan_variance=[compute_allan_variance( data,tau) for tau in tau_s]
plot((tau_s.*.5),(allan_variance),xlabel="log milliseconds",xaxis=:log, yaxis=:log,minorgrid=true, ylabel="allan variance (a.u)")
# plot((tau_s),(allan_variance))

# savefig(det_data_path[1:end-4]*"allan_variance")
# data= @view L[17000:40000]
# tau_s=1:floor(Int,length(data)/2)

# allan_variance=[compute_allan_variance( data,tau) for tau in tau_s]
# # tau_s=1:div(length(L_prime[4000:end]),2,RoundDown)    
# # allan_variance2=[compute_allan_variance_fast(L_prime[4000:end],tau) for tau in tau_s]
# plot(10*log10.(tau_s),10*log10.(allan_variance))
#  plot((tau_s),(allan_variance))
# f(x)=sin(x)
# x=1:0.001:1000
# data=f.(x)
# tau_s=1:floor(Int,length(data)/2)
# allan_variance=[compute_allan_variance( data,tau) for tau in tau_s]
# plot((tau_s),(allan_variance),xlabel="log milliseconds",xaxis=:log, yaxis=:log,minorgrid=true)
