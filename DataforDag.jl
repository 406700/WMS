## compare L and L_direct.

include("gas_experiment_functions.jl")
using LaTeXStrings, MAT
total_width_in_inches = 3.5
num_subplots = 1
subplot_width_in_inches = total_width_in_inches / num_subplots
subplot_height_in_inches = total_width_in_inches*4/5
fig, axs = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_height_in_inches))

#function to convert the temperature scale to wavenumber,
function map_range(x, in_min, in_max, out_min, out_max)
    return (x - in_min) / (in_max - in_min) * (out_max - out_min) + out_min
end

#function to calibrate the center wavelegnth
shift(x, shift) = x .+ shift

#function to 'calibrate' the wavenumber scale
function stretch(x_axis, scaling_factor)
    x_center = (maximum(x_axis) + minimum(x_axis)) / 2
    return x_center .+ scaling_factor .* (x_axis .- x_center)
end

hitran_path = "data/spectraplot/"
files=readdir(hitran_path)
hitran_path=hitran_path*files[1]
println("choosen hitran data: $hitran_path")

if true
    direct_path="data/20241211/direct_long.mat"
    mod_path = "data/20241211/20_long.mat" 
    # mod_path = "data/20250107/2khz_20.txt"

    # # Load data
    mod_L = load(mod_path[1:end-3] * "_L.jld2") 
    mod_x = load(mod_path[1:end-3] * "_xaxis.jld2") 

    direct_L = load(direct_path[1:end-3] * "sig_direct.jld2")
    direct_I0 = load(direct_path[1:end-3] * "ref_direct.jld2")
    direct_x = load(direct_path[1:end-3] * "_L_direct_temp.jld2")

    # Create the figure and axes
    fig, ax = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))

    for key in ["1"]#keys(direct_L)
        average_time=0.5 #seconds

        window_size = nearest_odd(average_time*1e6) 
        x,y = average_vector_in_chunks(direct_x[key],direct_L[key], window_size)
        _,I0 = average_vector_in_chunks(direct_x[key],direct_I0[key], window_size)
        y = y ./ I0

        window_size=nearest_odd(average_time*1e6/500)
        x2,y2=average_vector_in_chunks(mod_x[key],mod_L[key],window_size)
        normalization=1/maximum(y2)  


# Saving the data to a file called "figure_3_a_data.mat"
        data_to_save = Dict(
            "x_direct" => x,
            "y_direct_normalized" => y * normalization,
            "x_modulated" => x2,
            "y_modulated_normalized" => y2 * normalization
        )

        MAT.matwrite("figure_3_a_data.mat", data_to_save)

        data_to_save = Dict(
            "x_direct" => x,
            "y_direct" => y ,
            "x_modulated" => x2,
            "y_modulated" => y2 
        )

        MAT.matwrite("figure_3_a_data_not_normalized.mat", data_to_save)

    end
    

    ###############################################################################comparision of Lprime_over_L and -1/L( dα dν)
    using Dierckx
    #NB sign should correspong to derivative with respect to wavenumber not temp
    fig, ax = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_height_in_inches))
    p=Plots.plot()
    for average_time in [0.5]
        for key in ["1"]#keys(direct_L)
            # average_time=1 #seconds

            window_size = nearest_odd(average_time*1e6) 
            x,y = average_vector_in_chunks(direct_x[key],direct_L[key], window_size)
            _,I0 = average_vector_in_chunks(direct_x[key],direct_I0[key], window_size)
            y = y ./ I0

            #calculate the spline fit and its derivative of the normalized direct measurement
            spline=Spline1D(x, y, x[10:10:end-10];  k=3) #use knot vectors with larger spacing
            yprime=[derivative(spline,x_value) for x_value in x]
            L_prime=load(mod_path[1:end-3] * "_L_prime.jld2") 

            window_size=nearest_odd(average_time*2000)
            x2,y2=average_vector_in_chunks(mod_x[key],L_prime[key],window_size)

            # ax.plot(x,-yprime/maximum(-yprime),label="direct")
            # ax.plot(x2,y2/(maximum(y2)),label="wavelength modulated")
            data_to_save = Dict(
            "x_direct" => x,
            "yprime_direct_normalized" => -yprime / maximum(-yprime),
            "x_modulated" => x2,
            "yprime_modulated_normalized" => y2 / maximum(y2)
        )

        MAT.matwrite("figure_3_b_data_normalized.mat", data_to_save)
            GC.gc()

        end
    end
    Plots.savefig(p,"spie_figures/spline_vs_direct")

    # Finalize plot
    ax.legend()
    ax.set_xlabel(L"temp setpoint ($^\circ$ C)")  # Replace with appropriate label
    ax.set_ylabel(L"normalized $T^\prime$")  # Replace with appropriate label
    plt.tight_layout()
    plt.savefig("spie_figures/derivative_comparison.png",dpi=600)


    ####################################################################

    fig, ax = plt.subplots(2,1, figsize=(total_width_in_inches, subplot_height_in_inches))
    for average_time in collect(0.2:0.2:1)
        for key in ["1"]#keys(direct_L)
            # average_time=1 #seconds

            window_size = nearest_odd(average_time*2000) 
            x,y = average_vector_in_chunks(mod_x[key],mod_L[key], window_size)
            
            L_prime=load(mod_path[1:end-3] * "_L_prime.jld2") 
            x2,y2=average_vector_in_chunks(mod_x[key],L_prime[key],window_size)

            ax[1].plot(x,y)
            ax[2].plot(x2,y2)

            GC.gc()
        end
    end

    # Finalize plot
    # ax.legend()
    ax[1].set_xlabel(L"temp setpoint ($^\circ$ C)")  # Replace with appropriate label
    ax[1].set_ylabel(L"$T^\prime$")  # Replace with appropriate label
    plt.tight_layout()
    plt.savefig("spie_figures/transmittance_and_derivate.png",dpi=600)


    # spline2=Spline1D(rollx, -L_prime, rollx[200:400:end];  k=3)
    # Plots.plot(spline2.(rollx))



    # ax.legend()
    # ax.set_xlabel("temperature")  # Replace with appropriate label
    # ax.set_ylabel("scaled derivative, uncalibrated")  # Replace with appropriate label
    # plt.tight_layout()

    # plt.savefig("spie_figures/derivative_comparison.png",dpi=600)

    # x=1:1000
    # y=x.^2
    # sp1d=Spline1D(x, y;  k=3, bc="nearest", s=3.0)
    # Plots.plot(sp1d.(x))


    # Plots.plot(sp1d.(x[1:5000:end]))
    # Plots.plot!(y[1:5000:end])
    # 200/0.2


    # ###########################################################################hitran plot
# plot the numerical derivative of hitran data, converted to alpha


    det_data_paths=["data/20241211/20_long.mat"]
    # det_data_paths=["data/20250107/2khz_20_loss.mat"]


    key="1" #choose which scan to plot
    L=5.1 #length of the absorption cell
    # dwavenumber_to_dnu=2.99e10 #unit conversion
    Δνh=(-2.99e8*55e-12/2)/(1547e-9^2) # approximate value assumeing <h> is one and data for 20mv from first paper. factor of 2 from definition of chirp as half the total value

    #load hitran
    data, header = readdlm(hitran_path, ',', header = true) 
    hitran_x = data[:, 1]
    hitran_y = log.(data[:, 2])*-1/L #absorbance per cm from eq T=T0*exp(-alpha*L) => -log(T/T0)/L=alpha 
    hitran_xprime,hitran_yprime =numerical_derivative(hitran_x,hitran_y) #*100 to convert to d nu in meter.
    hitran_yprime=hitran_yprime.*1/2.99e10#convert to derivative with respect to nu

    #hitran frequency derivative. 
    

    xshift = 0.38
    xstretch = 0.8
    ystretch = 1.0

    #load precomputed experimental T^\prime/T
    L_prime_over_L=load(mod_path[1:end-3] * "Lprime_over_L.jld2") 
    average_time=0.5 #averaging time in seconds
    window_size = nearest_odd(average_time*2000) #calculate a window size equivalent to 0.5 seconds based on modulation rate of 2khz. 
    x,y=average_vector_in_chunks(mod_x[key],L_prime_over_L[key],window_size)
    y=-1/L*(y/-Δνh) #from eq 5 d_alpha/d_nu =-1/L T^\prime/T  deltanu<h>=T^\prime NB Δν_h=-1 in functions, so sign already accounted for but not magnitude. thus the - sign.

    x=map_range.(x,minimum(x),maximum(x),minimum(hitran_x),maximum(hitran_x)) #convert the x axis to inverse cm
    x = stretch(x, xstretch) #calibrate the axis 

    #scaling and plotting
    #NB inverse of the 
    scaling_factor=maximum(hitran_yprime)/maximum(-y)
    println("the 'extra' calibration factor related to error in Δν_h or hitran data is $scaling_factor")
    start_index=findfirst(x->x>=6462, hitran_xprime) #cut the hitran data to match the range of the experimental data

    fig, ax = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_height_in_inches)) #create plot
    # ax.plot(hitran_xprime[start_index:end], hitran_yprime[start_index:end], label = "HITRAN")
    # ax.plot(x.+xshift, -reverse(y)*scaling_factor,label=L"$T^\prime/T$ ")#, label = "Experimental Line $(i): $(det_data_path)")
    data_to_save = Dict(
            "x_wavenumber" => data[:, 1],
            "transmittance_51mm" => data[:, 1],
            "x_hertz" => hitran_xprime,
            "derivative_absorbance_per_cm_hz" => hitran_yprime,
            "absorbance_per_cm" => hitran_yprime
        )

        MAT.matwrite("Hitran.mat", data_to_save)

    # Finalize plot
    ax.legend()
    ax.set_xlabel(L"wavenumber $(cm^{-1})$")  # Replace with appropriate label
    ax.set_ylabel(L"$d\alpha/d\nu$ $(cm^{-1} Hz^{-1})$")  # Replace with appropriate label
    plt.tight_layout()
    plt.savefig("spie_figures/20241211_20_long_hitran_L_prime_over_L.png",dpi=600)

    ################################################################################################################################################3
    # with loss 


    fig, ax = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))
    fig2, ax2 = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))

    # ─────────────────────────────────────────────────────────────────────────────
    # Load your data

    mod_path    = "data/20241220/2khz_20.mat"
    ld_data_path = "data/20241220/2khz_20.txt"

    mod_L  = load(mod_path[1:end-3] * "_L.jld2")
    mod_x  = load(mod_path[1:end-3] * "_xaxis.jld2")
    L_prime_over_L = load("data/20241220/2khz_20_LprimeL.jld2")

    mod_path_loss = "data/20241220/2khz_20_loss.mat"
    ld_data_path_loss = "data/20241220/2khz_20_with_loss.txt"

    L_prime_over_L2 = load("data/20241220/2khz_LprimeL_loss.jld2")
    mod_L2 = load(mod_path_loss[1:end-3] * "_L.jld2")
    mod_x2 = load(mod_path_loss[1:end-3] * "_xaxis.jld2")

    L = 5.1
    # dwavenumber_to_dnu = 2.99e10
    # Δνh = -2.99e8*12e-9/(1547e-9^2)

    # ─────────────────────────────────────────────────────────────────────────────
    # Example loop (using key "1" only, as in your code)
    for key in ["1","2"]
        # Time-averaging setup
        average_time = 0.5  # seconds
        window_size = nearest_odd(average_time * 2000)

        # First figure data (dα/dν from L_prime_over_L)
        x,  y  = average_vector_in_chunks(mod_x[key],  L_prime_over_L[key],  window_size)
        x2_, y2 = average_vector_in_chunks(mod_x2[key], L_prime_over_L2[key], window_size)

        ax.plot(x,  y,  label="no loss")
        ax.plot(x2_, y2, label="~30% loss")

        # Second figure data (L from mod_L)
        x3, y3 = average_vector_in_chunks(mod_x[key],  mod_L[key],  window_size)
        x4, y4 = average_vector_in_chunks(mod_x2[key], mod_L2[key], window_size)

        ax2.plot(x3, y3, label="no loss")
        ax2.plot(x4, y4, label="30 percent loss")

        GC.gc()  # optional garbage collection
    end

    # ─────────────────────────────────────────────────────────────────────────────
    # Finalize and save first figure
    ax.legend()
    ax.set_xlabel(L"temp setpoint ($^\circ$ C)")  # or your preferred label
    ax.set_ylabel(L"$d\alpha/d\nu$") 
    fig.tight_layout()
    fig.savefig("spie_figures/loss_comparison_alpha.png", dpi=600)

    # ─────────────────────────────────────────────────────────────────────────────
    # Finalize and save second figure
    ax2.legend()
    ax2.set_xlabel(L"temp setpoint ($^\circ$ C)")  
    ax2.set_ylabel(L"$L$")  
    fig2.tight_layout()
    fig2.savefig("spie_figures/loss_comparison_T.png", dpi=600)


    ##########################################################################################@ 5 kHz


    fig, ax = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))
    fig2, ax2 = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))

    # ─────────────────────────────────────────────────────────────────────────────
    # Load your data

    mod_path    = "data/20241220/5khz_20.mat"
    ld_data_path = "data/20241220/5khz_20.txt"

    mod_L  = load(mod_path[1:end-3] * "_L.jld2")
    mod_x  = load(mod_path[1:end-3] * "_xaxis.jld2")
    L_prime_over_L = load("data/20241220/5khz_20_LprimeL.jld2")

    mod_path_loss = "data/20241220/5khz_20_loss.mat"
    ld_data_path_loss = "data/20241220/5khz_20_with_loss.txt"

    L_prime_over_L2 = load("data/20241220/5khz_LprimeL_loss.jld2")
    mod_L2 = load(mod_path_loss[1:end-3] * "_L.jld2")
    mod_x2 = load(mod_path_loss[1:end-3] * "_xaxis.jld2")

    L = 5.1
    # dwavenumber_to_dnu = 2.99e10
    # Δνh = -2.99e8*12e-9/(1547e-9^2)

    # ─────────────────────────────────────────────────────────────────────────────
    # Example loop (using key "1" only, as in your code)
    for key in ["1","2"]
        # Time-averaging setup
        average_time = 0.5  # seconds
        window_size = nearest_odd(average_time * 2000)

        # First figure data (dα/dν from L_prime_over_L)
        x,  y  = average_vector_in_chunks(mod_x[key],  L_prime_over_L[key],  window_size)
        x2_, y2 = average_vector_in_chunks(mod_x2[key], L_prime_over_L2[key], window_size)

        ax.plot(x,  y,  label="no loss")
        ax.plot(x2_, y2, label="~30% loss")

        # Second figure data (L from mod_L)
        x3, y3 = average_vector_in_chunks(mod_x[key],  mod_L[key],  window_size)
        x4, y4 = average_vector_in_chunks(mod_x2[key], mod_L2[key], window_size)

        ax2.plot(x3, y3, label="no loss")
        ax2.plot(x4, y4, label="30 percent loss")

        GC.gc()  # optional garbage collection
    end

    # ─────────────────────────────────────────────────────────────────────────────
    # Finalize and save first figure
    ax.legend()
    ax.set_xlabel(L"temp setpoint ($^\circ$ C)")  # or your preferred label
    ax.set_ylabel(L"$d\alpha/d\nu$") 
    fig.tight_layout()
    fig.savefig("spie_figures/loss_comparison_alpha_5khz.png", dpi=600)

    # ─────────────────────────────────────────────────────────────────────────────
    # Finalize and save second figure
    ax2.legend()
    ax2.set_xlabel(L"temp setpoint ($^\circ$ C)")  
    ax2.set_ylabel(L"$L$")  
    fig2.tight_layout()
    fig2.savefig("spie_figures/loss_comparison_T_5khz.png", dpi=600)


    ###########################################################################################20250107


    fig, ax = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))
    fig2, ax2 = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))

    # ─────────────────────────────────────────────────────────────────────────────
    # Load your data

    mod_path    = "data/20250107/2khz_20.mat"
    ld_data_path = "data/20250107/2khz_20.txt"

    mod_L  = load(mod_path[1:end-3] * "_L_prime.jld2")
    mod_x  = load(mod_path[1:end-3] * "_xaxis.jld2")
    L_prime_over_L = load("data/20250107/2khz_20.Lprime_over_L.jld2")

    mod_path_loss = "data/20250107/2khz_20_loss.mat"
    ld_data_path_loss = "data/20250107/2khz_20_loss.txt"

    L_prime_over_L2 = load("data/20250107/2khz_20_loss.Lprime_over_L.jld2")
    mod_L2 = load(mod_path_loss[1:end-3] * "_L_prime.jld2")
    mod_x2 = load(mod_path_loss[1:end-3] * "_xaxis.jld2")

    L = 5.1
    # dwavenumber_to_dnu = 2.99e10
    # Δνh = -2.99e8*12e-9/(1547e-9^2)

    # ─────────────────────────────────────────────────────────────────────────────
    # Example loop (using key "1" only, as in your code)
    for key in ["1","2"]
        # Time-averaging setup
        average_time = 0.5  # seconds
        window_size = nearest_odd(average_time * 2000)

        # First figure data (dα/dν from L_prime_over_L)
        x,  y  = average_vector_in_chunks(mod_x[key],  L_prime_over_L[key],  window_size)
        x2_, y2 = average_vector_in_chunks(mod_x2[key], L_prime_over_L2[key], window_size)

        ax.plot(x,  y,  label="no loss")
        ax.plot(x2_, y2, label="loss")

        # Second figure data (L from mod_L)
        x3, y3 = average_vector_in_chunks(mod_x[key],  mod_L[key],  window_size)
        x4, y4 = average_vector_in_chunks(mod_x2[key], mod_L2[key], window_size)

        ax2.plot(x3, y3, label="no loss")
        ax2.plot(x4, y4, label="loss")

        GC.gc()  # optional garbage collection
    end

    # ─────────────────────────────────────────────────────────────────────────────
    # Finalize and save first figure
    ax.legend()
    ax.set_xlabel(L"temp setpoint ($^\circ$ C)")  # or your preferred label
    ax.set_ylabel(L"$d\alpha/d\nu$") 
    fig.tight_layout()
    fig.savefig("spie_figures/loss_comparison_alpha_2025.png", dpi=600)

    # ─────────────────────────────────────────────────────────────────────────────
    # Finalize and save second figure
    ax2.legend()
    ax2.set_xlabel(L"temp setpoint ($^\circ$ C)")  
    ax2.set_ylabel(L"$L$")  
    fig2.tight_layout()
    fig2.savefig("spie_figures/loss_comparison_T_2025.png", dpi=600)
end

###########################################################################remaking conference plots with new data

direct_path="data/20250110/20250110_direct.mat"
mod_path = "data/20250109/2khz_20.mat" 
# mod_path = "data/20250107/2khz_20.txt"

# # Load data
mod_L = load(mod_path[1:end-3] * "_L.jld2") 
mod_x = load(mod_path[1:end-3] * "_xaxis.jld2") 

direct_L = load(direct_path[1:end-3] * "sig_direct.jld2")
direct_I0 = load(direct_path[1:end-3] * "ref_direct.jld2")
direct_x = load(direct_path[1:end-3] * "_L_direct_temp.jld2")

# Create the figure and axes
 fig, ax = plt.subplots(1, 1, figsize=(total_width_in_inches, subplot_height_in_inches))

 for key in ["1"]#keys(direct_L)
    average_time=0.5 #seconds

    window_size = nearest_odd(average_time*1e6) 
    x,y = average_vector_in_chunks(direct_x[key],direct_L[key], window_size)
    _,I0 = average_vector_in_chunks(direct_x[key],direct_I0[key], window_size)
    y = y ./ I0

    window_size=nearest_odd(average_time*1e6/500)
    x2,y2=average_vector_in_chunks(mod_x[key],mod_L[key],window_size)
    normalization=1/maximum(y2)  

    ax.plot(x, y*normalization, label="direct")
    ax.plot(x2,y2*normalization,label="wavelength modulated")

     GC.gc()
end

# Finalize plot
ax.legend()
ax.set_xlabel(L"temp setpoint ($^\circ$ C)")  # Replace with appropriate label
ax.set_ylabel(L"normalized $T$")  # Replace with appropriate label
plt.tight_layout()
plt.savefig("spie_figures/20mvmod_vs_direct_final.png",dpi=600)


using Dierckx
#NB sign should correspong to derivative with respect to wavenumber not temp
fig, ax = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_height_in_inches))
p=Plots.plot()
for average_time in [0.5]
    for key in ["1"]#keys(direct_L)
        # average_time=1 #seconds

        window_size = nearest_odd(average_time*1e6) 
        x,y = average_vector_in_chunks(direct_x[key],direct_L[key], window_size)
        _,I0 = average_vector_in_chunks(direct_x[key],direct_I0[key], window_size)
        y = y ./ I0

        #calculate the spline fit and its derivative of the normalized direct measurement
        spline=Spline1D(x, y, x[10:10:end-10];  k=3) #use knot vectors with larger spacing
        yprime=[derivative(spline,x_value) for x_value in x]
        L_prime=load(mod_path[1:end-3] * "_L_prime.jld2") 

        window_size=nearest_odd(average_time*2000)
        x2,y2=average_vector_in_chunks(mod_x[key],L_prime[key],window_size)

        ax.plot(x,-yprime/maximum(-yprime),label="direct")
        ax.plot(x2,y2/(maximum(y2)),label="wavelength modulated")

        GC.gc()

        Plots.plot!(p,x,y,label="direct")
        Plots.plot!(p,x,spline.(x),label="spline")

    end
end
Plots.savefig(p,"spie_figures/spline_vs_direct_final")

# Finalize plot
ax.legend()
ax.set_xlabel(L"temp setpoint ($^\circ$ C)")  # Replace with appropriate label
ax.set_ylabel(L"normalized $T^\prime$")  # Replace with appropriate label
plt.tight_layout()
plt.savefig("spie_figures/derivative_comparison_final.png",dpi=600)

#HITRAN
   hitran_path = "data/spectraplot/"
    files=readdir(hitran_path)
    hitran_path=hitran_path*files[1]
    println("choosen hitran data: $hitran_path")
    # det_data_paths=["data/20250107/2khz_20_loss.mat"]


    key="1" #choose which scan to plot
    L=5.1 #length of the absorption cell
    # dwavenumber_to_dnu=2.99e10 #unit conversion
    Δνh=(-2.99e8*55e-12/2)/(1547e-9^2) # approximate value assumeing <h> is one and data for 20mv from first paper. factor of 2 from definition of chirp as half the total value

    #load hitran
    data, header = readdlm(hitran_path, ',', header = true) 
    hitran_x = data[:, 1]
    hitran_y = log.(data[:, 2])*-1/L #absorbance per cm from eq T=T0*exp(-alpha*L) => -log(T/T0)/L=alpha 
    hitran_xprime,hitran_yprime =numerical_derivative(hitran_x,hitran_y) #*100 to convert to d nu in meter.
    hitran_yprime=hitran_yprime.*1/2.99e10#convert to derivative with respect to nu

    xshift = 0.38
    xstretch = 0.8
    ystretch = 1.0

    #load precomputed experimental T^\prime/T
    L_prime_over_L=load(mod_path[1:end-3] * "Lprime_over_L.jld2") 
    average_time=0.5 #averaging time in seconds
    window_size = nearest_odd(average_time*2000) #calculate a window size equivalent to 0.5 seconds based on modulation rate of 2khz. 
    x,y=average_vector_in_chunks(mod_x[key],L_prime_over_L[key],window_size)
    y=-1/L*(y/-Δνh) #from eq 5 d_alpha/d_nu =-1/L T^\prime/T  deltanu<h>=T^\prime NB Δν_h=-1 in functions, so sign already accounted for but not magnitude. thus the - sign.

    x=map_range.(x,minimum(x),maximum(x),minimum(hitran_x),maximum(hitran_x)) #convert the x axis to inverse cm
    x = stretch(x, xstretch) #calibrate the axis 

    #scaling and plotting
    #NB inverse of the 
    scaling_factor=maximum(hitran_yprime)/maximum(-y)
    println("the 'extra' calibration factor related to error in Δν_h or hitran data is $scaling_factor")
    start_index=findfirst(x->x>=6462, hitran_xprime) #cut the hitran data to match the range of the experimental data

    fig, ax = plt.subplots(1,1, figsize=(total_width_in_inches, subplot_height_in_inches)) #create plot
    ax.plot(hitran_xprime[start_index:end], hitran_yprime[start_index:end], label = "HITRAN")
    ax.plot(x.+xshift, -reverse(y)*scaling_factor,label=L"$T^\prime/T$ ")#, label = "Experimental Line $(i): $(det_data_path)")

    # Finalize plot
    ax.legend()
    ax.set_xlabel(L"wavenumber $(cm^{-1})$")  # Replace with appropriate label
    ax.set_ylabel(L"$d\alpha/d\nu$ $(cm^{-1} Hz^{-1})$")  # Replace with appropriate label
    plt.tight_layout()
    plt.savefig("spie_figures/20250109_L_prime_over_L_final.png",dpi=600)


