function analytic_n0(x)
	return exp(x)-1
end

function analytic_n1(x)
	return -x-1+exp(x)
end

function analytic_n2(x)
	return -.5*x^2-x-1+exp(x)
end

function analytic_n3(x)
	return -.5*x^2-x-1+exp(x)-x^3/6.0
end

println("Analytic solution")
println(analytic_n0(10^(-4)))
println(analytic_n1(10^(-4)))
println(analytic_n2(10^(-4)))
println(analytic_n3(10^(-4)))

function integrand0(u,x)
	return u^0*exp(0)
end

function integrand1(u,x)
	return u*exp(-u*x)
end

function integrand2(u,x)
	return u^2*exp(-u*x)
end

function integrand3(u,x)
	return u^3*exp(-u*x)
end

function front_term(x,n)
	return x^(n+1)*exp(x)/factorial(n)
end

function newton_two(f, a, b, N, x)
	delta = (b-a)/(1.0*N)
	h = delta/2.0
	sum = h*f(a, x)/3+h*f(b, x)/3
	for n in 1:N-1
		sum+=(h/3)*(2*f(a+2*n*h, x)+4*f(a+(2*n-1)*h, x))
	end
	sum+= 4*h*f(a+(2*N-1)*h, x)/3
	return sum
end

x = 10^(-4)

println("Simpsons estimates")

#Evaluate integral for n=0
integral_term = newton_two(integrand0, 0, 1, 1000, x)
res0 = front_term(x, 0)*integral_term
println(res0)

#Evaluate integral for n=1
integral_term = newton_two(integrand1, 0, 1, 1000, x)
res1 = front_term(x, 1)*integral_term
println(res1)


#Evaluate integral for n=2
integral_term = newton_two(integrand2, 0, 1, 1000, x)
res2 = front_term(x, 2)*integral_term
println(res2)

#Evaluate integral for n=3
integral_term = newton_two(integrand3, 0, 1, 1000, x)
res3 = front_term(x, 3)*integral_term
println(res3)


function relative_error(numerical, true_sol)
	numerator = abs(numerical-true_sol)
	return numerator/(1.0*true_sol)
end

println("Relative Errors")

#For analytical results:
i_0_exact = 1.0000500016667083e-4
x0=relative_error(analytic_n0(10^(-4)), i_0_exact)
y0=relative_error(res0, i_0_exact)

i_1_exact = 5.0001666708334167e-9
x1=relative_error(analytic_n1(10^(-4)), i_1_exact)
y1=relative_error(res1, i_1_exact)

i_2_exact = 1.6667083341666806e-13
x2=relative_error(analytic_n2(10^(-4)), i_2_exact)
y2=relative_error(res2, i_2_exact)

i_3_exact = 4.1667500013877985e-18
x3=relative_error(analytic_n3(10^(-4)), i_3_exact)
y3=relative_error(res3, i_3_exact)


x = [x0,x1,x2,x3]
y = [y0,y1,y2,y3]
n=[0,1,2,3]

using PyPlot
figure()
semilogy(n, x , label=" Analytical Solution")
semilogy(n, y,  label="Simpson's Rule")
legend()
xlabel("N")
ylabel("Relative Error")
title("Problem 3 - Relative Error for integrals I_n")
savefig("problem3_relerr.png")
print("done")







