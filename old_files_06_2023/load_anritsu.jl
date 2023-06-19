using DelimitedFiles
using Plots

cd("/home/m/OneDrive/Experimental_Data/Chirped_laser_fall/")
cd("20221114_15mv/")

OSA_files=readdir("./OSA/")
osc_files=readdir("./oscilloscope/")

readdlm("./oscilloscope//"*osc_files[1],',',skipstart=11)

function read_osc_data(files)
    osc_header=readdlm("./oscilloscope/"*files[1],',')[1:11,:]
    osc_data=Array{Any}(undef,length(files))
    for i in 1:length(files)
        osc_data[i]=readdlm("./oscilloscope/"*files[i],',',skipstart=11)
    end
    return osc_data, osc_header
end

function read_osa_data(files)
    osa_header=readdlm("./OSA/"*files[1],',')[1:27,:]
    osa_data=Array{Any}(undef,length(files))
    for i in 1:length(files)
        osa_data[i]=readdlm("./OSA/"*files[i],',',skipstart=29)
    end
    return osa_data, osa_header
end


foo,_=read_osa_data(OSA_files)

foo,_=read_osc_data(osc_files)


#fit a lorentzian to the spectra and find the center wavelength
function OSA_center_wavelength()
end

function treat_osc_data(data,pulse_length)
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>=0,time_axis)
    data_array=data[t_0:pulse_length,3]
    time_axis=time_axis[t_0:pulse_length]

    for i in 2:length(files)
        data=readdlm(path*files[i],',',skipstart=12)
        t_0=findfirst(x->x>0,data[:,1])
        data_array=hcat(data_array,data[t_0:pulse_length,3])
    end

    return time_axis,data_array, Temps
end
