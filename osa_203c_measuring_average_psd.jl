
using DelimitedFiles
using Plots
plotlyjs()
#20231219 commenting out below 2 lines
#cd("/home/m/OneDrive/Experimental_Data/2023_04_26_measuring_chirp_on_osa")
#cd("/home/m/OneDrive/Experimental_Data/20230503_measuring_chirp_on_osa/round_7/")
cd("/home/m/OneDrive/Experimental_Data/WMS_final_paper/")

files=readdir()

f=files[11]
data=readdlm(f,',',skipstart=2)
datastart=findfirst(x->x=="[Data]",data)[1]
header=data[1:datastart,:]
data=data[(datastart+1):(end-1),:]
y_unit=findfirst(x->x=="#YAxisUnit",header)[1]
y_unit=header[y_unit,2]

start=1546.8
stop=1547.005
start_index=argmin(abs.(data.-start))[1]
stop_index=argmin(abs.(data.-stop))[1]
data=data[start_index:stop_index,:]
plot(data[:,1],data[:,2],xlabel="wavelength (nm)",ylabel=y_unit)
#findfirst(x->x<data[1,2],data[10:end,2])

freqxpower=0
for i in 1:length(data[1:end,1])
    global freqxpower=freqxpower+data[i,1]*data[i,2]
end
average=freqxpower/sum(data[:,2])
vline!([average])

#savefig("average_without_tail.png")