using Statistics
using Plots
using DelimitedFiles
using LsqFit
using EasyFit
using Interpolations

plotlyjs()
include("S_measurement_data.jl")
#default(show=true)

#interactively choose pulse length. Return the data and the temperature of the TEC for each measurement.
function read_data(path)
    files=readdir(path)
    Temps=get_temps(files) #gets the temperature of the TEC controller on which the grating is mounted from the file name.
    #data_array=Array{Float64,2}(undef,length(files))

    #read the first file and set the length of time axis from user input
    data=readdlm(path*files[1],',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis)
    plt=plot()
    gui(plot(plt,data[t_0:end,3]))
    pulse_length=0

    while true
        println("enter pulse length")
        pulse_length=readline()
        pulse_length=parse(Int,pulse_length)
        plot!(plt,data[1:pulse_length,3],show=true)
        println("accept value enter y or n")
        input=readline()
        if lowercase(input) == "y"
            println("processing...")
            break
        end
    end

    data_array=data[t_0:pulse_length,3]
    time_axis=time_axis[t_0:pulse_length]

    for i in 2:length(files)
        data=readdlm(path*files[i],',',skipstart=12)
        t_0=findfirst(x->x>0,data[:,1])
        data_array=hcat(data_array,data[t_0:pulse_length,3])
    end

    return time_axis,data_array, Temps
end

#accepts pulse length as argument, instead of running interactively
function read_data(path,pulse_length)
    files=readdir(path)
    Temps=get_temps(files)

    #data_array=Array{Float64,2}(undef,length(files))

    #read the first file and set the length of time axis from user input
    data=readdlm(path*files[1],',',skipstart=12)#time,channel1, channel,2volts seconds
    time_axis=data[:,1]*1e6 #convert to us
    t_0=findfirst(x->x>0,time_axis)
    data_array=data[t_0:pulse_length,3]
    time_axis=time_axis[t_0:pulse_length]

    for i in 2:length(files)
        data=readdlm(path*files[i],',',skipstart=12)
        t_0=findfirst(x->x>0,data[:,1])
        data_array=hcat(data_array,data[t_0:pulse_length,3])
    end

    return time_axis,data_array, Temps
end

#read the temperature of the TEC from the file path
function get_temps(files)
    Temps = zeros(length(files))
    for i in 1:length(files)
        Temps[i] = parse(Float64, files[i][1:4]) #take the first 3 digits of file name and convert to float64
    end
    return Temps
end


function calculate_S(t,data, t1, t2,T)
    S = copy(T)
    t1=findfirst(x->abs(x-t1)<.00001,t) #NB could lead to offset. Could also just enter index instead of time.
    t2=findfirst(x->abs(x-t2)<.00001,t)
    println("t1=$t1, t2 = $t2")
    for i = 1:length(S)
        p_minus = data[t1,i]
        p_plus = data[t2,i]
        S[i] = (p_plus - p_minus) / (p_plus + p_minus)
    end
    return S
end

#interactively choose two time values and return S
function plot_S_interactive(t,data,T)
    gui(plot(t,data,label=transpose(T)))
    t1=0
    t2=0
    while true
        println("enter t1")
        t1=readline()
        t1=parse(Float64,t1)
        println("enter t2")
        t2=readline()
        t2=parse(Float64,t2)

        println("accept values y or n?")
        y_n=readline()
        if lowercase(y_n) == "y"
            println("calculating S")
            break
        end
    end
    println("t1=$t1, t2=$t2")
    S=calculate_S(t,data,t1,t2,T)
    gui(plot(T,S))
    return t1,t2,S
end


experiment_list=ChirpedLaser_S_Data.s_experiment_list
t,f_t,T=read_data(experiment_list[2].path)
t1,t2,S=plot_S_interactive(t,f_t,T)

plot(t,f_t,label=transpose(T))

s=calculate_S(t,f_t,50,60,T)
plot(T,s)
plot(s,show=true)
