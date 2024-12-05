include("gas_experiment_functions.jl") 
load_det_data=false

global const sampling_rate=1e6
global const mod_rate=2e3
global const data_points_per_period=sampling_rate/mod_rate
global const data_points_per_half_period=round(Int,data_points_per_period/2)

det_data_path="data/20241202/03bar_stab.mat" #gas briefly disconnected and hose removed, before restarting the experiment with only a few minutes flow time. Assumed gas may still be present.

if load_det_data == false

    ref_chan,sig_chan,time_chan=read_det_data(det_data_path,0.995) #trigger level, determines when laser powers on and off
    GC.gc()
   
    ref_chan_dict = Dict{Int, Vector{Float64}}()
    sig_chan_dict = Dict{Int, Vector{Float64}}()
    time_chan_dict = Dict{Int, Vector{Float64}}()
    number_scans=1
    GC.gc()

    for index in 1:number_scans
        key = index  # Using integers as keys
        ref_chan_dict[key] = ref_chan
        sig_chan_dict[key] = sig_chan
        time_chan_dict[key] = time_chan

    end
    ref_chan=nothing
    sig_chan=nothing
    GC.gc()
 # Save all six dictionaries to disk
 
    open(det_data_path[1:end-3]*"bin", "w") do io
        serialize(io, (
            ref_chan_dict,
            sig_chan_dict,
            time_chan_dict
           
        ))
    end
else
    # Load all six dictionaries from disk
    (
        ref_chan_dict,
        sig_chan_dict,
        time_chan_dict
    ) = open(det_data_path[1:end-3]*"bin", "r") do io
        deserialize(io)
    end
end


L_dict = Dict{Int, Array{Float64}}()
L_prime_dict = Dict{Int, Array{Float64}}()
trig_dict = Dict{Int, Tuple{Float64, Float64}}()

# Loop over all keys in the dictionaries
function get_normalization_coefficient(ref_chan,sig_chan)
    

    downsample = 100
    # Find reference trigger level
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)

    # Find trigger indices
    trig_indices = find_pulse_trig_points(ref_chan, ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

    # Trim channels
    ref_chan, sig_chan, trig_indices = trim_channels(ref_chan, sig_chan, trig_indices, 2)
    return mean(ref_chan[trig_indices[1]:trig_indices[3]]),mean(sig_chan[trig_indices[1]:trig_indices[3]])
end

ref_norm,sig_norm=get_normalization_coefficient(ref_chan_dict[1],sig_chan_dict[1])

for i in keys(ref_chan_dict)
    println("key=$i")
    # Extract the channels for the current key
    ref_chan = ref_chan_dict[i]
    sig_chan = sig_chan_dict[i]
    time_chan = time_chan_dict[i]
    downsample = 100
    # Find reference trigger level
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)

    # Find trigger indices
    trig_indices = find_pulse_trig_points(ref_chan, ref_trig_level_intercept, ref_trig_level_slope, data_points_per_half_period)

    # Trim channels
    ref_chan, sig_chan, trig_indices = trim_channels(ref_chan, sig_chan, trig_indices, 2)

    # Normalize channels
    # ref_chan, sig_chan = normalize_channels(ref_chan, sig_chan, trig_indices)
    ref_chan=ref_chan/ref_norm
    sig_chan=sig_chan/sig_norm

    # Recalculate the reference trigger level for trimmed and normalized data
    ref_trig_level_slope, ref_trig_level_intercept = find_ref_trigger_level(ref_chan, downsample)
    trig_dict[i]= (ref_trig_level_slope, ref_trig_level_intercept)
    ref_chan_dict[i]= ref_chan 
    sig_chan_dict[i] = sig_chan
    time_chan_dict[i] = time_chan
    # Run garbage collection
    GC.gc()

    # Calculate L and L_prime for the current scan
    L_dict[i],L_prime_dict[i],_,_ = calculate_L_L_prime_for_scan(ref_chan, sig_chan, trig_indices, ref_trig_level_slope, ref_trig_level_intercept)

    # Store results in the dictionary
    
end
for key in keys(L_dict)
    if iseven(key)
        L_dict[key]=reverse(L_dict[key])  
        L_prime_dict[key]=reverse(L_prime_dict[key])  
    end
end

p=plot()
for key in keys(L_prime_dict)
    plot!(p,L_prime_dict[key])
    # plot!(p,L_dict[key])
end
display(p)
savefig(det_data_path[1:end-3]*"_L_prime.png")

p=plot()
for key in keys(L_dict)
    plot!(p,L_dict[key])
    # plot!(p,L_dict[key])
end
display(p)

savefig(det_data_path[1:end-3]*"_L.png")

save(det_data_path[1:end-3]*"_L.jld2",Dict(string(key) => value for (key, value) in L_dict))
save(det_data_path[1:end-3]*"_L_prime.jld2",Dict(string(key) => value for (key, value) in L_prime_dict))


