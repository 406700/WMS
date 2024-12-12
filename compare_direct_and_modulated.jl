
using JLD2,Plots,RollingFunctions

direct_path="data/20241211/direct_long.mat"
direct_path="data/20241211/direct_long.txt"

mod_path="data/20241211/20_long.mat" 
mod_path="data/20241211/20_long.txt"

mod_L=load(mod_path[1:end-3] * "_L.jld2") 
mod_x= load(mod_path[1:end-3] * "_xaxis.jld2") 

direct_L=load(direct_path[1:end-3]*"_L_direct.jld2")
direct_x=load(direct_path[1:end-3]*"_L_direct_temp.jld2")

p=plot()

for key in ["2"]#keys(direct_L)[1]
    window_size=5001
    y=rollmean(direct_L[key],window_size)
    half_w = (window_size - 1) ÷ 2
    rollx = direct_x[key][half_w+1:end-half_w]    
    plot!(p,rollx[1:100:end],y[1:100:end],label="direct")

    window_size=9
    y=rollmean(direct_L[key],window_size)
    half_w = (window_size - 1) ÷ 2
    rollx = mod_x[key][half_w+1:end-half_w]    
    plot!(p,rollx,y,label="mod")

end
display(p)