include("gas_experiment_functions.jl")
using LaTeXStrings

total_width_in_inches = 3.5
num_subplots = 3
subplot_width_in_inches = total_width_in_inches / num_subplots
subplot_height_in_inches = total_width_in_inches * 4 / 5

fig, axs = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))
fig2, axs2 = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))

# ────────────────────────────────────────────────────────────────────────
# Load your data

mod_path    = "data/20250109/2khz_20.mat"
ld_data_path = "data/20250109/2khz_20.txt"

mod_L  = load(mod_path[1:end-3] * "_L.jld2")
mod_L_prime  = load(mod_path[1:end-3] * "_L_prime.jld2")

mod_x  = load(mod_path[1:end-3] * "_xaxis.jld2")
# L_prime_over_L = load("data/20250107/2khz_20.Lprime_over_L.jld2")

mod_path_loss = "data/20250109/2khz_20_loss.mat"
ld_data_path_loss = "data/20250109/2khz_20_loss.txt"

# L_prime_over_L2 = load("data/20250107/2khz_20_loss.Lprime_over_L.jld2")
mod_L2 = load(mod_path_loss[1:end-3] * "_L.jld2")
mod_L_prime2 = load(mod_path_loss[1:end-3] * "_L_prime.jld2")
mod_x2 = load(mod_path_loss[1:end-3] * "_xaxis.jld2")


mod_path_loss_t = "data/20250109/2khz_20_loss_t.mat"
ld_data_path_loss_t = "data/20250109/2khz_20_loss_t.txt"

# L_prime_over_L2 = load("data/20250107/2khz_20_loss.Lprime_over_L.jld2")
mod_L3 = load(mod_path_loss_t[1:end-3] * "_L.jld2")
mod_L_prime3 = load(mod_path_loss_t[1:end-3] * "_L_prime.jld2")
mod_x3 = load(mod_path_loss[1:end-3] * "_xaxis.jld2")

L = 5.1

# ────────────────────────────────────────────────────────────────────────


# fig, axs = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))
# fig2, axs2 = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))

# for key in ["1", "2"]
#     # Time-averaging setup
#     average_time = 0.5  # seconds
#     window_size = nearest_odd(average_time * 2000)

#     # First figure data (dα/dν from L_prime_over_L)
#     x, y = average_vector_in_chunks(mod_x[key], L_prime_over_L[key], window_size)
#     axs[1].plot(x, y, label="")

#     x2, y2 = average_vector_in_chunks(mod_x2[key], L_prime_over_L2[key], window_size)
#     axs[1].plot(x2, y2, label="loss")

#     x3, y3 = average_vector_in_chunks(mod_x_new[key], L_prime_over_L_new[key], window_size)
#     axs[2].plot(x3, y3, label="")

#     # Second figure data (L from mod_L)
#     x4, y4 = average_vector_in_chunks(mod_x_new2[key], L_prime_over_L_new2[key], window_size)
#     axs[2].plot(x4, y4, label="loss")

#     x5, y5 = average_vector_in_chunks(mod_x[key], mod_L[key], window_size)
#     axs2[1].plot(x5, y5, label="")

#     x5, y5 = average_vector_in_chunks(mod_x2[key], mod_L2[key], window_size)
#     axs2[1].plot(x5, y5, label="loss")

#     x6, y6 = average_vector_in_chunks(mod_x_new[key], mod_L_new[key], window_size)
#     axs2[2].plot(x6, y6, label="")

#     x6, y6 = average_vector_in_chunks(mod_x_new2[key], mod_L_new2[key], window_size)
#     axs2[2].plot(x6, y6, label="")
#     GC.gc()  # optional garbage collection
# end

# # ────────────────────────────────────────────────────────────────────────
# # Finalize and save first figure
# axs[1].set_xlabel(L"temp setpoint ($^\circ$ C)")
# axs[1].set_ylabel(L"$d\alpha/d\nu$")
# for i in 1:2
#     axs[i].legend()
#     axs[i].set_title("Dataset $(i)")
# end

# fig.tight_layout()
# fig.savefig("spie_figures/loss_comparison_alpha_2025.png", dpi=600)

# # Finalize and save second figure
# axs2[1].set_xlabel(L"temp setpoint ($^\circ$ C)")
# axs2[1].set_ylabel(L"$L$")
# for i in 1:2
#     axs2[i].legend()
#     axs2[i].set_title("Dataset $(i)")
# end

# fig2.tight_layout()
# fig2.savefig("spie_figures/loss_comparison_T_2025.png", dpi=600)

#----------------------------------------------------------------------------


fig, axs = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))
fig2, axs2 = plt.subplots(1, 2, figsize=(total_width_in_inches, subplot_height_in_inches))

for key in ["1", "2"]
    # Time-averaging setup
    average_time = 0.5  # seconds
    window_size = nearest_odd(average_time * 2000)

    # First figure data (dα/dν from L_prime_over_L)
    x, y = average_vector_in_chunks(mod_x[key], mod_L_prime[key]./mod_L[key], window_size)
    axs[1].plot(x, y, label="")

    x2, y2 = average_vector_in_chunks(mod_x2[key], mod_L_prime2[key]./mod_L2[key], window_size)
    axs[1].plot(x2, y2, label="loss")


    x, y = average_vector_in_chunks(mod_x[key], mod_L_prime[key]./mod_L[key], window_size)
    axs[2].plot(x, y, label="")

    x3, y3 = average_vector_in_chunks(mod_x3[key], mod_L_prime3[key]./mod_L3[key], window_size)
    axs[2].plot(x3, y3, label="")

    # Second figure data (L from mod_L)
   
    x5, y5 = average_vector_in_chunks(mod_x[key], mod_L[key], window_size)
    axs2[1].plot(x5, y5, label="")

    x5, y5 = average_vector_in_chunks(mod_x2[key], mod_L2[key], window_size)
    axs2[1].plot(x5, y5, label="loss")

    x6, y6 = average_vector_in_chunks(mod_x[key], mod_L[key], window_size)
    axs2[2].plot(x6, y6, label="")

    x6, y6 = average_vector_in_chunks(mod_x3[key], mod_L3[key], window_size)
    axs2[2].plot(x6, y6, label="")
    GC.gc()  # optional garbage collection
end

# ────────────────────────────────────────────────────────────────────────
# Finalize and save first figure
axs[1].set_xlabel(L"temp setpoint ($^\circ$ C)")
axs[1].set_ylabel(L"$d\alpha/d\nu$")
for i in 1:2
    axs[i].legend()
    axs[i].set_title("Dataset $(i)")
end

fig.tight_layout()
fig.savefig("spie_figures/loss_comparison_alpha_202509.png", dpi=600)

# Finalize and save second figure
axs2[1].set_xlabel(L"temp setpoint ($^\circ$ C)")
axs2[1].set_ylabel(L"$L$")
for i in 1:2
    axs2[i].legend()
    axs2[i].set_title("Dataset $(i)")
end

fig2.tight_layout()
fig2.savefig("spie_figures/loss_comparison_T_202509.png", dpi=600)
