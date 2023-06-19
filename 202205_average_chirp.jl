
using DelimitedFiles
using Plots
plotlyjs()

function center_spectrum(data)
    center_wavelength=1546.898 #wavelength meter round 5
    half_span=40
    idx=findfirst(x->x>center_wavelength,data[:,1])
    data=data[(idx-half_span):(idx+half_span),:]
end

function voltage_to_laser_current(modulation_amplitude)
    m=0.150
    ILD_set=26.6#22+16*m
    modulation_amplitude=modulation_amplitude/2
    ILD = ILD_set - modulation_amplitude * m
end

voltage_to_laser_current(32)
cd("/home/m/OneDrive/Experimental_Data/20230503_measuring_chirp_on_osa/round_7/")
files=readdir()[1:end]


#for i in 1:length(files)

    f=files[2]
    data=readdlm(f,',',skipstart=2)
    datastart=findfirst(x->x=="[Data]",data)[1]
    header=data[1:datastart,:]
    data=data[(datastart+1):(end-1),:]
    data=center_spectrum(data)
    y_unit=findfirst(x->x=="#YAxisUnit",header)[1]
    y_unit=header[y_unit,2]

    #data=data[1:15,:]
    #plot!(data[1:end,1],data[1:end,2],xlabel="wavelength (nm)",ylabel=y_unit)
    #findfirst(x->x<data[1,2],data[10:end,2])

    freqxpower=0
    for i in 1:length(data[1:end,1])
        global freqxpower=freqxpower+data[i,1]*data[i,2]
    end
    average=freqxpower/sum(data[:,2])
    #vline!([average])
    #averages[i]=average
#end
#savefig("average_without_tail.png")
p=plot()
titles=[1,3,5,8,11,14,17,20,23,26,29,32]
averages=zeros(length(files))
for i in 1:length(files)

    f=files[i]
    data=readdlm(f,',',skipstart=2)
    datastart=findfirst(x->x=="[Data]",data)[1]
    header=data[1:datastart,:]
    data=data[(datastart+1):(end-1),:]
    data=center_spectrum(data)
    y_unit=findfirst(x->x=="#YAxisUnit",header)[1]
    y_unit=header[y_unit,2]

    #data=data[1:15,:]
    plot!(p,data[1:end,1],data[1:end,2],xlabel="wavelength (nm)",ylabel=y_unit)#,title=titles[i])
    #findfirst(x->x<data[1,2],data[10:end,2])
    freqxpower=0
    for i in 1:length(data[1:end,1])
       freqxpower=freqxpower+data[i,1]*data[i,2]
    end
    average=freqxpower/sum(data[:,2])
    vline!([average])
    
    averages[i]=average
end
gui(p)


