using LinearAlgebra
function gradient_hardcode(xk, Q, b)
	top = xk[1]*Q[1]+xk[2]*Q[3]-b[1]
	bot = xk[1]*Q[2]+xk[2]*Q[4]-b[2]
	return [top, bot]
end

function mag(vec1)
	return (vec1[1]^2+vec1[2]^2)^(.5)
end

function update(xk, alpha, Q, b)
	g = gradient_hardcode(xk, Q, b)
	return (xk[1]-alpha*g[1], xk[2]-alpha*g[2])
end

function make_alpha(xk, Q, b)
	g = gradient_hardcode(xk, Q, b)
	numerator = g[1]^2+g[2]^2
	s_top = g[1]*Q[1]+g[2]*Q[3]
	s_bot = g[1]*Q[2]+g[2]*Q[4]
	denominator = s_top*g[1]+s_bot*g[2]
	return numerator/denominator
end

function iterate_steepest_descent(x0, a_func, b, Q)
	xk = x0
	alpha = a_func(xk,Q, b)
	next_x = update(xk, alpha, Q, b)
	rel_change = mag([next_x[1]-xk[1], next_x[2]-xk[2]])
	other_change = 10^(-8.0)*mag(next_x)
	allx = [(x0[1],x0[2])]
	while rel_change > (other_change)
		push!(allx, next_x)
		res=evaluate_f(next_x,Q,b)
		xk = next_x
		alpha = a_func(xk,Q, b)
		next_x = update(xk, alpha, Q, b)
		rel_change = mag([next_x[1]-xk[1], next_x[2]-xk[2]])
		other_change = 10^(-8.0)*mag(next_x)
	end
	return allx
end

function rel_err(allx, fixedpt)
	final = []
	for i in 1:size(allx)[1]
		numerator = mag((allx[i][1]-fixedpt[1],allx[i][2]-fixedpt[2]))
		total = numerator/norm(fixedpt)
		push!(final, total)
	end
	return final
end

function a_constant(xk,b,Q)
	return .1
end

function evaluate_f(xk,Q,b)
	quad_term = xk[1]^2*Q[1]+xk[2]*xk[1]*Q[3]+xk[2]*xk[1]*Q[2]+xk[2]^2*Q[4]
	b_term = b[1]*xk[1]+b[2]*xk[2]
	return quad_term-b_term
end

b = [2.0 3.0]
Q = [4.0 -2.0; -2.0 3.0]
x0 = [-4.0 -3.0]
xs = iterate_steepest_descent(x0, a_constant, b, Q)

xfinal = xs[size(xs)[1]]
println(xfinal)
res=evaluate_f(xfinal,Q,b)
println(res)


using PyPlot

# Example data
n = 100
x = LinRange(-6, 6, n)
y = LinRange(-6, 6, n)

xgrid = repeat(x', n, 1)
ygrid = repeat(y, 1, n)
zgrid = xgrid .* ygrid

z = rand(n, n)
for i = 1:n
    for j = 1:n
    	z(i,j) = evaluate_f((xgrid[i],ygrid[j]),Q,b)
	end
end
# This next line is not strictly necessary, but being able to access the figure's handle could be useful to set parameters.
fig = figure()  

# Get the current axis (GCA) and store it as a variable
ax = gca()

# Create the contour plot (cp)
cp = ax[:contour](xgrid, ygrid, zgrid, colors="black", linewidth=2.0)

# Label the contour plot
ax[:clabel](cp, inline=1, fontsize=10)

x=[]
y=[]
println(size(xs))
for x1 in 1:size(xs)[1]
	push!(x,xs[x1][1])
	push!(y,xs[x1][2])
end
plot(x,y)

title("Contour Plot")
savefig("problem5_constantalpha.png")
println("done")

figure()
out_err = rel_err(xs, (1.5,2.0))
x = LinRange(0, size(out_err)[1], size(out_err)[1])
semilogy(x, out_err)
xlabel("Number of Iterations")
ylabel("Relative Error")
title("Constant Alpha Steepest Descent Error")
savefig("problem5_constalpha_err.png")


#Optimal Alpha Code
xs = iterate_steepest_descent(x0, make_alpha, b, Q)
xfinal = xs[size(xs)[1]]
println(xfinal)
res=evaluate_f(xfinal,Q,b)
println(res)

println(size(xs))

using PyPlot

# Example data
n = 100
x = LinRange(-6, 6, n)
y = LinRange(-6, 6, n)

xgrid = repeat(x', n, 1)
ygrid = repeat(y, 1, n)
zgrid = xgrid .* ygrid

z = rand(n, n)
for i = 1:n
    for j = 1:n
    	z(i,j) = evaluate_f((xgrid[i],ygrid[j]),Q,b)
	end
end
# This next line is not strictly necessary, but being able to access the figure's handle could be useful to set parameters.
fig = figure()  

# Get the current axis (GCA) and store it as a variable
ax = gca()

# Create the contour plot (cp)
cp = ax[:contour](xgrid, ygrid, zgrid, colors="black", linewidth=2.0)

# Label the contour plot
ax[:clabel](cp, inline=1, fontsize=10)

x=[]
y=[]
println(size(xs))
for x1 in 1:size(xs)[1]
	push!(x,xs[x1][1])
	push!(y,xs[x1][2])
end
plot(x,y)

title("Contour Plot")
savefig("problem5_optimalalpha.png")
println("done")


figure()
out_err = rel_err(xs, (1.5,2.0))
x = LinRange(0, size(out_err)[1], size(out_err)[1])
semilogy(x, out_err)
xlabel("Number of Iterations")
ylabel("Relative Error")
title("Optimal Alpha Steepest Descent Error")
savefig("problem5_optalpha_err.png")






using LinearAlgebra
function newton_nd(f, J, x0; N=1)
    x = x0
    println(x)
    for i=1:N
    	r=(J(x) \ f(x))
        x = [x[1]-r[1], x[2]-r[2]]
    end 
    println(x)
    return x
end

function get_tolerance(newtonND, f, J, x0)
    n=1
    x_N=newtonND(f,J,x0; N=1)
    x_list = []
    while mag(f(x_N))>10^(-8.0)
    	n+=1
    	push!(x_list, x_N)
    	x_N=newtonND(f,J,x0;N=n)
    end
    return x_N, x_list, n
end

function rel_err(allx, fixedpt)
	final = []
	for i in 1:size(allx)[1]
		numerator = mag((allx[i][1]-fixedpt[1],allx[i][2]-fixedpt[2]))
		total = numerator/norm(fixedpt)
		push!(final, total)
	end
	return final
end



f(xk) = [xk[1]*Q[1]+xk[2]*Q[3]-b[1]; xk[1]*Q[2]+xk[2]*Q[4]-b[2]]
J(x) = [
    (Q[1])   		(Q[3])
    (Q[2])			(Q[4])
]

#x0=[-10,-10]
x_N, x_list, iters  = get_tolerance(newton_nd, f, J, x0)
#println(evaluate_f([x_N[1],x_N[2]],Q,b))
#println(x_list)
#println(iters)
out_err = rel_err(x_list, [x_N[1],x_N[2]])
x=[]
println(size(x_list))
println(out_err)
for i in 1:iters
	push!(x,i)
end
print(x)
#semilogy(x, out_err)
#title("Relative Error")
#savefig("problem5_newton_err.png")



