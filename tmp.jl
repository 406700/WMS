using Interpolations,FiniteDiff,Plots
# Define data points
x = collect(range(-10,stop=10,length=1000))
y = x.^3

# # Create a linear interpolation
# interpolated_line = LinearInterpolation(x, y, extrapolation_bc=NaN)

# # Define the function using the interpolation
# f(x_val) = interpolated_line(x_val)

# # Use finite difference to compute the derivative at all original data points
# line_derivative =[ FiniteDiff.finite_difference_derivative(f, x_point) for x_point in x] 
line_derivative=zeros(1000)
for i in 2:length(x)-1
    # Central difference formula
    line_derivative[i]  = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    # println("Derivative at x = $(x[i]) using finite difference: $line_derivative")
end

plot(x,y)
plot!(x,line_derivative)