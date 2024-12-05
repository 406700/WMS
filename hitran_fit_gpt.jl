using DelimitedFiles, Interpolations, Plots, RollingFunctions, JLD2,FiniteDiff,Statistics

# Utility Functions
shift(x, shift) = x .+ shift
function numerical_derivative(x,y)

   line_derivative = zeros(Float64, length(y))  
   for i in 2:Int(length(x)-1)
        # Central difference formula
        line_derivative[i]  = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    end
    return x,line_derivative
end
function stretch(x_axis, scaling_factor)
    x_center = (maximum(x_axis) + minimum(x_axis)) / 2
    return x_center .+ scaling_factor .* (x_axis .- x_center)
end

# Function 1: Compare with HITRAN Data
function compare_hitran(hitran_path, det_data_paths, direct_measurement)
    files = readdir(hitran_path)
    data, header = readdlm(hitran_path * files[1], ',', header = true)
    x = data[:, 1]
    y = data[:, 2]

    p = plot(x, y, label = "HITRAN Data")
    χ_shift = 0.8
    plot!(p, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN χ $χ_shift")

    for (i, det_data_path) in enumerate(det_data_paths)
        exp_line = load(det_data_path[1:end-3] * "_L.jld2")["1"] 
        if direct_measurement == true
            exp_line = rollmean(exp_line, 100)[1:500:end]
        end

        exp_x = collect(range(x[1], stop = x[end], length = length(exp_line)))
        # xshift = 0.1 
        # xstretch = 0.7
        # ystretch = 1.0
        # yshift = 0#0.01
        xshift = 0#0.1 
        xstretch = 1#0.7
        ystretch = 1.0
        yshift = 0#0.01

        exp_x = stretch(exp_x, xstretch)
        exp_x = shift.(exp_x, xshift)
        exp_line2 = stretch(exp_line, ystretch)
        exp_line2 = shift.(exp_line, yshift)

        plot!(p, exp_x, exp_line2)#, label = "Experimental Line $(i): $(det_data_path)")
    end

    return p
end
function compare_hitran2(hitran_path, det_data_paths, direct_measurement)
    files = readdir(hitran_path)
    data, header = readdlm(hitran_path * files[4], ',', header = true)
    x = data[:, 1]
    y = reverse(data[:, 2])

    p = plot(x, y, label = "HITRAN Data")
    χ_shift = 0.8
    plot!(p, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN χ $χ_shift")

    for (i, det_data_path) in enumerate(det_data_paths)
        exp_line = load(det_data_path[1:end-3] * "_L.jld2")["1"] 
        if direct_measurement == true
            exp_line = rollmean(exp_line, 100)[1:500:end]
        end

        exp_x = collect(range(x[1], stop = x[end], length = length(exp_line)))
        # xshift = 0.1 
        # xstretch = 0.7
        # ystretch = 1.0
        # yshift = 0#0.01
        xshift = -0.1#0.1 
        xstretch = 0.9
        ystretch = 1.0
        yshift = 0.02#+0.01#0.01

        exp_x = stretch(exp_x, xstretch)
        exp_x = shift.(exp_x, xshift)
        exp_line2 = stretch(exp_line, ystretch)
        exp_line2 = shift.(exp_line, yshift)

        plot!(p, exp_x, exp_line2)#, label = "Experimental Line $(i): $(det_data_path)")
    end

    return p
end
# Function 2: Compare Direct and Derived Lines
function compare_direct_and_derived_lines(direct_path, mod_paths, pres_paths)
    key = "1"
    direct_line = load(direct_path[1:end-3] * "_L.jld2")[key]
    direct_line = rollmean(direct_line, 100)[1:500:end]  # Smooth and downsample

    p = plot()
    x_direct = collect(range(start = 1, stop = length(direct_line), length = length(direct_line)))

    # Plot modulation data
    for mod in keys(mod_paths)
        mod_line = load(mod_paths[mod][1:end-3] * "_L.jld2")[key]
        println(size(mod_line))
        global x_direct = collect(range(start = 1, stop = length(mod_line), length = length(direct_line)))
        plot!(p, 1:length(mod_line), mod_line, label = "Modulation $(mod)")
    end

    # Plot direct line
    plot!(p, x_direct, direct_line, label = "Direct Measurement")

    # Plot pressure data
    for (pressure, path) in pres_paths
        pres_line = load(path[1:end-3] * "_L.jld2")[key]
        plot!(p, 1:length(pres_line), pres_line, label = "Pressure $(pressure) bar",legend=false)
    end
    
    return p
end

# Function 1: Compare with HITRAN Data
function compare_hitran_derivative(hitran_path, det_data_paths, direct_measurement)
    data, header = readdlm(hitran_path, ',', header = true)
    x = data[:, 1]
    y = data[:, 2]
    p=plot()
    p2 = plot()
    χ_shift = 1

    plot!(p2, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN χ $χ_shift")
    x,y =numerical_derivative(x,y)

    plot!(p, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN χ $χ_shift")
    max_deriv=maximum(y)
    xshift = 0.4
    xstretch = 0.8
    ystretch = 1.0
    yshift = 0#0.015
    for (i, det_data_path) in enumerate(det_data_paths)
        exp_line = load(det_data_path[1:end-3] * "_L_prime.jld2")["1"] 
        # if direct_measurement == true
        #     exp_line = rollmean(exp_line, 100)[1:500:end]
        # end

        exp_x = collect(range(x[1], stop = x[end], length = length(exp_line)))
        
        exp_x = stretch(exp_x, xstretch)
        exp_x = shift.(exp_x, xshift)
        exp_line2 = stretch(exp_line, ystretch)
        # exp_line2 = shift.(exp_line, yshift)
        # exp_line2=-exp_line2#nb just reversed the sign to match
        plot!(p, exp_x, reverse(exp_line2)/maximum(exp_line2)*max_deriv)#, label = "Experimental Line $(i): $(det_data_path)")
    end
    #plot the line
    for (i, det_data_path) in enumerate(det_data_paths)
        exp_line = load(det_data_path[1:end-3] * "_L.jld2")["1"] 
        # if direct_measurement == true
        #     exp_line = rollmean(exp_line, 100)[1:500:end]
        # end

        exp_x = collect(range(x[1], stop = x[end], length = length(exp_line)))
        exp_x = stretch(exp_x, xstretch)
        exp_x = shift.(exp_x, xshift)
        exp_line2 = stretch(exp_line, ystretch)
        exp_line2 = shift.(exp_line, yshift)

        plot!(p2, exp_x, reverse(exp_line2))#nb_reverse, label = "Experimental Line $(i): $(det_data_path)")
    end

    return p,p2
end
function compare_hitran_derivative_with_adjusted_scale(hitran_path, det_data_paths)
    data, header = readdlm(hitran_path, ',', header = true)
    x = data[:, 1]
    y = data[:, 2]
    p=plot(xlabel="wavenumber",ylabel="uncalibrated derivative (a.u.)")
    p2 = plot(xlabel="wavenumber",ylabel="Transmittance")
    χ_shift = 1

    plot!(p2, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN χ $χ_shift")
    x,y =numerical_derivative(x,y)

    plot!(p, x, (1 .-((1 .-y)* χ_shift)), label = "HITRAN χ $χ_shift")
    max_deriv=maximum(y)
    xshift = 0.4
    xstretch = 0.38
    ystretch = 1.0
    yshift = 0#0.015

    for (i, det_data_path) in enumerate(det_data_paths)
        exp_line = load(det_data_path[1:end-3] * "_L_prime.jld2")["3"] 
        exp_x= load(det_data_path[1:end-3] * "_xaxis.jld2")["3"] 
        exp_x=map_range.(exp_x,minimum(exp_x),maximum(exp_x),minimum(x),maximum(x))
        # if direct_measurement == true
        #     exp_line = rollmean(exp_line, 100)[1:500:end]
        # end
      
        exp_x = stretch(exp_x, xstretch)
        # exp_x = shift.(exp_x, xshift)
        exp_line2 = stretch(exp_line, ystretch)
        # exp_line2 = shift.(exp_line, yshift)
        # exp_line2=-exp_line2#nb just reversed the sign to match
        plot!(p, exp_x.+xshift, reverse(exp_line2)/maximum(exp_line2)*max_deriv)#, label = "Experimental Line $(i): $(det_data_path)")
   
        # if direct_measurement == true
        #     exp_line = rollmean(exp_line, 100)[1:500:end]
        # end
        exp_line2 = load(det_data_path[1:end-3] * "_L.jld2")["3"] 
        # exp_line2 = stretch(exp_line, ystretch)
        # exp_line2 = shift.(exp_line, yshift)
        plot!(p2, exp_x.+xshift, reverse(exp_line2))
    end

    return p,p2
end
function map_range(x, in_min, in_max, out_min, out_max)
    return (x - in_min) / (in_max - in_min) * (out_max - out_min) + out_min
end

# Main Code
hitran_path = "data/spectraplot/"
files=readdir(hitran_path)
hitran_path=hitran_path*files[1]
# det_data_paths = ["data/20241128/gas_no_mod.mat"]  
# det_data_paths=["data/20241128/gas_mod10.mat"]# ,"data/20241128/gas_mod10.mat"]
# det_data_paths=["data/20241202/01bar_long.mat"]# ,"data/20241128/gas_mod10.mat"]
det_data_paths=["data/20241202/03bar_scan.mat"]

direct_measurement = false
# p1 = compare_hitran(hitran_path, det_data_paths , direct_measurement)
# p1 = compare_hitran2(hitran_path, det_data_paths , direct_measurement)
# p1,p2=compare_hitran_derivative(hitran_path, det_data_paths,false)
p1,p2=compare_hitran_derivative_with_adjusted_scale(hitran_path, det_data_paths)

savefig(p1,det_data_paths[1][1:end-3]*"hitran_L_prime.png") #nb for vector of paths
savefig(p2,det_data_paths[1][1:end-3]*"hitran_L.png")

display(p2)

## Direct measurement paths
# mod_paths = Dict(
#     30 => "data/20241125/10min_mod2.mat",
#     20 => "data/20241128/gas_mod20.mat",
#     10 => "data/20241128/gas_mod10.mat"
# )

# pres_paths = Dict(
#     1 => "data/20241129/01bar.mat",
#     2 => "data/20241129/02bar.mat"
# )

# direct_path = "data/20241128/gas_no_mod.mat"
# p2 = compare_direct_and_derived_lines(direct_path, mod_paths, pres_paths)
# display(p2)
