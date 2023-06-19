using DelimitedFiles
using Statistics
using Plots
using LsqFit
using EasyFit

L_prime=readdlm("/home/m/OneDrive/Experimental_Data/20230612_stability/L_prime_corrected")
foo=zeros(length(L_prime))
for i in 1:999
    foo[i]=L_prime[i]
end

##################################################################
#2nd order
L_prime=foo*1e12
t=[range(0,step=.001,length=999);]
p0=[1.0,1.0,1.0,1.0] #x0,gamma,amplitude
@. model(x,p)=p[1]+p[2]*(x-p[4])+p[3]*(x-p[4])^2
fit=curve_fit(model,t,L_prime,p0)
@. model(x)=fit.param[1]+fit.param[2]*(x-fit.param[4])+fit.param[3]*(x-fit.param[4])^2



plot(t,fit.resid*1e-12)
plot(t,model.(t),label="fit",xlabel="time [S]",ylabel="L_prime*10^12")
plot!(t,L_prime,label="data",title="measurement stability data vs fit")
ylims!(-4.2,-3.7)
mean_error=std(fit.resid*1e-12)

########################################################3
#linear 
# linear_fit =fitlinear(t,L_prime)
# plot(linear_fit.x,linear_fit.y)
# plot!(t,L_prime)
# mean_error=std(linear_fit.residues)*1e-12

##error in derivative/(slope of derivative function with lambda)
L_prime_with_mod=readdlm("L_prime_with_fp_modulation")*1e12
t=[t;]
slope=fitlinear(t[210:220],L_prime_with_mod[210:220])
plot(t,L_prime_with_mod)
plot!(slope.x,slope.y)
slope=slope.a*1e-12

measurement_error=mean_error/(slope*0.5/720e-12)*1e12 #pm

