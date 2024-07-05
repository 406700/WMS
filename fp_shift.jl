using Plots
plotlyjs()
using EasyFit
v=[2.7,2.8,2.9,3,3.1,3.6,4.1,4.6,5.1]
λ=[1544.604,1545.0143,1545.412,1545.800,1546.184,1548.042,1549.873,1551.59,1553.3]
plot(v,λ)
scatter!(v,λ)

fit=fitlinear(v,λ)
plot!(fit.x,fit.y)
fit.a*0.3
λ[4]-λ[1]
(1.09-1.19)/1.19
1.09#nm per volt, but slope is not completely linear, thus this should be different on different days. I can just deal with the error here?

fwhm=mean([216,217,219,220,221,210,219,219]) #don't know how accurate my measurements where, since I felt it was closer to 220pm.
fit=fitlinear(v[1:5],λ[1:5])
fit.a