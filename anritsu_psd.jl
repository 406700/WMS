
using DelimitedFiles
using Plots
plotlyjs()

###########26ma
dir1=("data/chirp_calibration/current_26_T15_mod_70_10/")
dir2=("data/chirp_calibration/current_26_T185_mod_70_10/")
dir3=("data/chirp_calibration/current_26_T21_mod_70_10/")

files1=readdir(dir1)
files2=readdir(dir2)
files3=readdir(dir3)

#############23ma
dir1=("data/chirp_calibration/current_23_T21_mod_70_0/")
# dir2=("data/chirp_calibration/current_26_T185_mod_70_10/")
# dir3=("data/chirp_calibration/current_26_T21_mod_70_10/")

files1=readdir(dir1)
files2=readdir(dir2)
files3=readdir(dir3)

##########

dir1=("data/chirp_calibration/current_285_T_mod_70_0/")
# dir2=("data/chirp_calibration/current_26_T185_mod_70_10/")
# dir3=("data/chirp_calibration/current_26_T21_mod_70_10/")

files1=readdir(dir1)
files2=readdir(dir2)
files3=readdir(dir3)

###########################################################
dir1=("data/wavelength_calibration_285/")
files1=readdir(dir1)

colors = [:red, :blue, :green, :orange, :purple,:yellow,:red,:blue,]
averagepower=zeros(length(files1))
p=plot()

for (i,f) in enumerate(files1) 
    data=readdlm(dir1*f,skipstart=28,',')
    plot!(p,data[:,1],data[:,2])#,color=colors[i])
    datalength=length(data[:,1])
    freqxpower=zeros(datalength)
    for i in 1:datalength
        freqxpower[i]=(data[i,1])*data[i,2]
    end
    power=sum(data[:,2])
    average=sum(freqxpower)/power
    averagepower[i]=average
    # vline!([average])#,color=colors[i])
end
display(p)
plot(avg,averagepower.-averagepower[end])

# # delta_lamba=collect(10:10:100)*1e-9
# # delta_nu=[2.99e8.*x/1547e-9^2 for x in delta_lamba]

# # plot(delta_nu)

