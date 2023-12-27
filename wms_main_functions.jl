
function L0(i_p,i_m,Δim)

    1/(2*(1-Δim^2))*(i_p-Δim*i_m)
end

function L0_prime(i_p,i_m,Δν_h,Δim)
   1/(Δν_h*2*(1-Δim^2))*(i_m-Δim*i_p)
end

function average_values(FP,indices,sample_periods,averages,T,ts_range,I_0)
    i_plus=zeros(sample_periods,averages)
    i_minus=zeros(sample_periods,averages)
    FP_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        FP_array[:,i]=FP[indices[i]:Int(indices[i]+T-1)]
    end
    for i in 1:averages
        ts=ts_range[i]
        i_plus[:,i]=FP_array[ts,:].+FP_array[Int(ts+T/2),:]
        i_minus[:,i]=FP_array[ts,:].-FP_array[Int(ts+T/2),:]
        # sample_plus=range(ts,step=T,length=sample_periods) #defines evenly spaced samples. Need to correct for drift in the trigger.
        # sample_minus=sample_plus.+Int(T/2)
        # i_plus[:,i]=FP[sample_plus].+FP[sample_minus]
        # i_minus[:,i]=FP[sample_plus].-FP[sample_minus]
    end
    i_plus=mean(i_plus,dims=2)
    i_minus=mean(i_minus,dims=2)
    return i_plus/I_0,i_minus/I_0
end

#treat the direct LD the same as the modulated data, except average for each period. 
function LD_average_values(L_direct,indices,sample_periods,T)
    LD_array=zeros(sample_periods)
    for i in 1:sample_periods
        LD_array[i]=mean(L_direct[indices[i]:Int(indices[i]+T-1)])
    end
   return LD_array
end

function dλ_to_dν(dI_dλ,x)
    λ= x.*1e-9 
    dI_dν=zeros(length(dI_dλ))

    dλ_dν=(λ.^2)./c
    ν=c./(λ)
    #c/λ^2
    for i in 1:length(dλ_dν)
        dI_dν[i]=dI_dλ[i]*1e-9*dλ_dν[i] #1e-9 to convert the derivative in units of nm to units of meters. 
    end
    return dI_dν, ν
end

#this really just prepares the data by formatting it the same way and dividing into each pulse.
function chirp_estimate(FP,LD,L_direct,sample_periods,T,indices)
    #chirp_h = zeros(sample_periods)
    #chirp_l = zeros(sample_periods)
    FP_array=zeros(T,sample_periods)
    direct_array=zeros(T,sample_periods)
    LD_array=zeros(T,sample_periods)
    for i in 1:sample_periods
        FP_array[:,i]=FP[indices[i]:Int(indices[i]+T-1)]./LD[indices[i]:Int(indices[i]+T-1)]
        LD_array[:,i]=LD[indices[i]:Int(indices[i]+T-1)]
        direct_array[:,i]=L_direct[indices[i]:Int(indices[i]+T-1)]
    end
    return FP_array,direct_array,LD_array
    
end

function get_chirp_value(chirp, x,L_direct_avg,λ)
    # Check if x is within the bounds of the array 'L_direct_avg'
    if x < 1 || x > length(L_direct_avg)
        error("Starting point x is out of bounds in vector L_direct_avg")
    end

    # Find max and min values and their positions in vector 'chirp'
     max_value, max_position = findmax(chirp)
     min_value, min_position = findmin(chirp)

    # Determine which one (max or min) occurs first in 'chirp'
    if max_position < min_position
        y, z = max_value, min_value
    else
        y, z = min_value, max_value
    end

    #y=mean(chirp[30000:31000])
    #z=mean(chirp[61250:62250])
    # Initialize indices for y and z in vector 'L_direct_avg'
    index_of_y = 1
    index_of_z = 1
    # Search backward for y in 'L_direct_avg'
    for i in x:-1:1
        if abs(L_direct_avg[i] - y) <= .0035 #calculated from spacing in L_direct_avg.
            index_of_y = i
            break
        end
    end

    # Search forward for z in 'L_direct_avg'
    for i in x:length(L_direct_avg)
        if abs(L_direct_avg[i] - z) <= .0035
            index_of_z = i
            break
        end
    end

    return abs(λ[index_of_y]-λ[index_of_z])
end

function get_chirp_estimate(FP_array, direct_array, LD_array,L_direct_avg,λ, N)
    chirp=zeros(N)
    for i in 1:N
        chirp[i]=get_chirp_value(FP_array[:,i],i,L_direct_avg,λ)
    end
    return chirp
end

function lorentzian_fit(λ,L_direct_avg)
        ################################################################ lorentzian fit 
    @. lor_model(x,p)= p[3]*(1/2*p[2])/(pi*(x-p[1])^2+(1/2*p[2])^2)#lorentzian model
    @. lor_derivative_model(x,p)= p[3]*(-16(x-p[1])*p[2])/(pi*4(x-p[1])^2+p[2]^2)^2#lorentzian model
    x=λ*1e9
    y=L_direct_avg
    p0=[1547,0.2,1.0] #initial guess
    lor_fit=curve_fit(lor_model,x,y,p0)

    L_direct_fit=lor_model(x,lor_fit.param)# 
    L_prime_direct=lor_derivative_model(x,lor_fit.param) #dI/dlambda nm
    return L_direct_fit, L_prime_direct
end

function convert_to_dv(L_prime_direct,λ)
    #convert to derivative to dI/dν from dI/d_lambda
    #direct_derivative, ν = dλ_to_dν(L_direct_fit,x) #fit params and lambda in nm 

    dI_dν=zeros(length(L_prime_direct))
    dλ_dν=(λ.^2)./c
    ν=c./(λ)
        #c/λ^2
    for i in 1:length(dλ_dν)
            dI_dν[i]=L_prime_direct[i]*1e9*dλ_dν[i] #1e-9 to convert the derivative in units of nm to units of meters. 
    end
    return dI_dν
end

function get_fwhm(λ, L,L_direct_avg)
        ################################################################# 
    #calculate the fwhm
    ind1=findfirst(x->x>(maximum(L)/2),L)
    ind2=findfirst(x->x==maximum(L),L)
    fwhm_method=abs(2*(λ[ind1]-λ[ind2]))
    ind1=findfirst(x->x>(maximum(L_direct_avg)/2),L_direct_avg)
    ind2=findfirst(x->x==maximum(L_direct_avg),L_direct_avg)
    fwhm_direct=abs(2*(λ[ind1]-λ[ind2]))
    return fwhm_method, fwhm_direct
end



function trigger_values(LD,T,sample_periods,search_width=1000)
    LD_trigger=LD.-mean(LD)
    #search_width=1000 #based on visually checking first and last LD pulse
    step_size=Int(T/2-search_width/2)
    indices=zeros(sample_periods*2)
    search_start=1
    for i in 1:length(indices) 
        next_ind=findmin(abs.(LD_trigger[search_start:search_start+step_size]))[2]#find the minimum element in the expected range 
        indices[i]=Int(search_start+next_ind)
        search_start=search_start+next_ind+step_size
    end
    indices=Int.(indices[1:2:(2*sample_periods-1)])
    return indices
end
function offset_trigger_values(LD,T,sample_periods,search_width=1000)
    LD_trigger=LD.-mean(LD)
    #search_width=1000 #based on visually checking first and last LD pulse
    step_size=Int(T/2-search_width/2)
    indices=zeros(sample_periods*2)
    search_start=1
    for i in 1:length(indices) 
        next_ind=findmin(abs.(LD_trigger[search_start:search_start+step_size]))[2]#find the minimum element in the expected range 
        indices[i]=Int(search_start+next_ind)
        search_start=search_start+next_ind+step_size
    end
    #indices=Int.(indices[1:2:(2*sample_periods-1)]) removed.
    return Int.(indices)
end

function find_rising_trigger(LD, indices, window_size)
    for index in indices
        ##check the slope around the point
        if LD[index-window_size]<LD[index+window_size]
            return index
        end    
    end
    return -1
end



