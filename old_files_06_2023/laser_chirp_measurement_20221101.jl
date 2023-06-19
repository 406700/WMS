using DelimitedFiles
using Plots
using LsqFit
using EasyFit
plotlyjs()
##NB oscilloscope spectrum analyzer channels are hard coded, so the indices would need to be changed if things are hooked up differently.

function read_osc_data(files)
    osc_header=readdlm("./oscilloscope/"*files[1],',')[1:12,:]
    osc_data=Array{Any}(undef,length(files))
    for i in 1:length(files)
        osc_data[i]=readdlm("./oscilloscope/"*files[i],',',skipstart=12)
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

#fit a lorentzian to the spectra and find the center wavelength

#cut the time to a single pulse
function time_osc_data(data,pulse_length)
    time_axis=data[:,1]
    pulse_index=findfirst(x->x>=pulse_length,time_axis)
    t_0=findfirst(x->x>=0,time_axis)
    time_axis=time_axis[t_0:pulse_index]
    data=hcat(time_axis,data[t_0:pulse_index,2:end])
end

#return the parameters of a lorentzian fit to the OSA data
function fit_osa_lorentzian(data)
        # p= amplitude,width,center
        model(x,p)=p[1]*p[2] ./((p[2]^2) .+(x .-p[3]).^2)
        p0=[145e-6,0.3,1547.6]
        xdata=data[:,1]
        ydata = data[:,2]
        fit=curve_fit(model,xdata,ydata,p0)
        return fit.param

     osa_fit_params
end


function loss_calibration(Osc)
    #This function normalizes the photodetector power measurements. # Assuming the transmission is equivalent to 1, this gives the differences in the losses in the system. 
    #losses= (p(t)*fp_max_transmission*connector_loss)/(p(t)*connector_loss_2).
    #Takes the laser transmission through the FP and divides it by the max #p(t)/p_0(t) for the time t where the laser wavelength is exactly at the peak of the transmission spectrum. 
    P_max=maximum(Osc[:,2])#maximum of sensor trace
    max_index=findfirst(x->x>=P_max,Osc[:,2])
    loss=Osc[max_index,2]/Osc[max_index,4] #p(t)/p_0(t)
end

function chirp(OSC,OSA,losses)
    #OSC[:,2] = laser->fp->det 
    #OSC[:,4] = laer->det (p_0(T))
    # takes the photodetector voltages(OSC_), the FP transmission spectrum(OSA), and the losses (see description in loss function)
    normalized_transmission=OSC[:,2]./(OSC[:,4].*losses) #fp_transmission/p_0(t)*losses. 
    #normalized_transmission=movingaverage(normalized_transmission,1000).x
    max_transmission=findmax(normalized_transmission)[2]#index of the peak of the lorentzian fp spectrum.
    #normalize the transmission spectrum on the OSA of the FP to be 1 at the maximum.
    normalize_constant=findmax(OSA[:,2])[1]
    max_index=findmax(OSA[:,2])[2] #peak of the spectrum
    FP_transmission=OSA[:,2]./normalize_constant

    #convert the p(t) to \lambda(t) 
    #find where the power transmitted throught the fp is equivalent to the fp transmission spectrum.
    #function is cut in two sections to simplify the determination of which side of the transmission peak it is on.
    #NB This assumes the transmission starts to the left of the peak and increases in wavelength. Won't work for second hald of cycle.

    function difference_1(x)
        ind=findmin(abs.(x.-FP_transmission[1:max_index]))[2]
        OSA[1:max_index][ind,1]
    end
    function difference_2(x)
        ind=findmin(abs.(x.-FP_transmission[max_index:end]))[2]
        OSA[max_index:end][ind,1]
    end
    λ_of_t=difference_1.(normalized_transmission[1:max_transmission-1])
    λ_of_t_2=difference_2.(normalized_transmission[max_transmission:end])
    return vcat(λ_of_t,λ_of_t_2)
end

function moving_average(OSC,n)
    column_2=movingaverage(OSC[:,2],n).x
    column_4=movingaverage(OSC[:,4],n).x
    return hcat(OSC[:,1],column_2,OSC[:,3],column_4)
end


function initialize()
    @time OSA_data,_=read_osa_data(OSA_files)
    # @time osa_fit_params=fit_osa_lorentzian(OSA_data)
    osa_fit_params=Array{Any}(undef,length(OSA_data))
    for i in 1:length(osa_fit_params)
        osa_fit_params[i]=fit_osa_lorentzian(OSA_data[i])
    end
    OSC_data,_=read_osc_data(osc_files)
    return OSA_data,OSC_data,osa_fit_params
end

cd("/home/m/OneDrive/Experimental_Data/Chirped_laser_fall/")
cd("20221123_45mv/")

OSA_files=readdir("./OSA/")
osc_files=readdir("./oscilloscope/")

OSA_data,OSC_data,osa_fit_params=initialize()


#calculate the chirp
selected_OSC_data=time_osc_data(OSC_data[8],500e-6)
selected_OSC_data=moving_average(selected_OSC_data,200)
losses=loss_calibration(selected_OSC_data)
λ_chirp=chirp(selected_OSC_data,OSA_data[8],losses)
#plot(λ_chirp)

#OSA_data[8][348,1]
#plot(OSC_data[7][:,1],OSC_data[7][:,4])
#plot!(OSC_data[8][:,1],OSC_data[8][:,2])

#plot 
#a=plot()
#for i in OSC_data
#    plot!(a,i[:,1],i[:,2:4])
#end
#gui(a)

#foo=time_osc_data(OSC_data[1][1],1000e-6)
#plot(foo[:,1],foo[:,2])
#plot(OSA_data[8][:,1],OSA_data[8][:,2]./0.0001561)


# function calculate_S(t,data, t1, t2,T)
#     S = copy(T)
#     t1=findfirst(x->abs(x-t1)<.00001,t) #NB could lead to offset. Could also just enter index instead of time.
#     t2=findfirst(x->abs(x-t2)<.00001,t)
#     println("t1=$t1, t2 = $t2")
#     for i = 1:length(S)
#         p_minus = data[t1,i]
#         p_plus = data[t2,i]
#         S[i] = (p_plus - p_minus) / (p_plus + p_minus)
#     end
#     return S
# end

# #interactively choose two time values and return S
# function plot_S_interactive(t,data,T)
#     gui(plot(t,data,label=transpose(T)))
#     t1=0
#     t2=0
#     while true
#         println("enter t1")
#         t1=readline()
#         t1=parse(Float64,t1)
#         println("enter t2")
#         t2=readline()
#         t2=parse(Float64,t2)

#         println("accept values y or n?")
#         y_n=readline()
#         if lowercase(y_n) == "y"
#             println("calculating S")
#             break
#         end
#     end
#     println("t1=$t1, t2=$t2")
#     S=calculate_S(t,data,t1,t2,T)
#     gui(plot(T,S))
#     return t1,t2,S
# end
# P_max=maximum(OSA[:,2])
# λ_max_index=findfirst(x->x>=P_max,OSA[:,2])
# λ_max=OSA[λ_max_index,2
