
using JLD2,Plots,RollingFunctions


# direct_path="data/20241211/dierct_13.mat"
direct_path="data/20241211/direct_long"
# direct_path="data/20241211/no_gas_long.mat"
# # direct_path="data/20241211/no_gas_long.txt"
mod_path="data/20241211/20_long.mat" 
mod_path="data/20241211/20_long.txt"


mod_L=load(mod_path[1:end-3] * "_L.jld2") 
mod_x= load(mod_path[1:end-3] * "_xaxis.jld2") 

direct_L=load(direct_path[1:end-3]*"sig_direct.jld2")
direct_I0=load(direct_path[1:end-3]*"ref_direct.jld2")

direct_x=load(direct_path[1:end-3]*"_L_direct_temp.jld2")

p=plot()

for key in keys(direct_L)
    window_size=5001
    y=rollmean(direct_L[key],window_size)
    I0=rollmean(direct_I0[key],window_size)
    half_w = (window_size - 1) ÷ 2
    rollx = direct_x[key][half_w+1:end-half_w]  
    x= rollx[1:100:end]
    y=y./I0
    y=y[1:100:end]
    plot!(p,x,y,label="direct")

    # window_size=9
    # y=rollmean(mod_L[key],window_size)
    # half_w = (window_size - 1) ÷ 2
    # rollx = mod_x[key][half_w+1:end-half_w]    
    # plot!(p,rollx,y,label="mod")
    # GC.gc()
end
display(p)
# savefig(p,mod_path[1:end-3]*"_corrected_scans_full.png")

