function newton_zero(f,a,b, N)
	delta = (b-a)/(1.0*N)
	sum=0.0
	for i in 1:N-1
		sum+= delta*f(a+i*delta)
	end
	return sum
end 

function newton_one(f,a,b, N)
	delta = (b-a)/(1.0*N)
	sum=.5*delta*(f(a)+f(b))
	for i in 1:N-1
		sum+= delta*f(a+i*delta)
	end
	return sum
end

function newton_two(f, a, b, N)
	delta = (b-a)/(1.0*N)
	h = delta/2.0
	sum = h*f(a)/3+h*f(b)/3
	for n in 1:N-1
		sum+=(h/3)*(2*f(a+2*n*h)+4*f(a+(2*n-1)*h))
	end
	sum+= 4*h*f(a+(2*N-1)*h)/3
	return sum
end

x = [10,100,1000,10000,100000,1000000,10000000]

function error(psum,exact_soln)
	num = abs(psum-exact_soln)
	return num/exact_soln
end

I_a_exact=2.5193079820307612557

function I_a(x)
	return exp(cos((x+1)^2+2*sin(4*x+1)))
end

y_0 = [error(newton_zero(I_a,0,pi,n), I_a_exact) for n in x]
y_1 = [error(newton_one(I_a,0,pi,n), I_a_exact) for n in x]
y_2 = [error(newton_two(I_a,0,pi,n), I_a_exact) for n in x]

using PyPlot
figure()
loglog(x,y_0,label=L"0-th order")
loglog(x,y_1,label= L"1-st order")
loglog(x,y_2, label = L"2-nd order")
legend()
xlabel("N")
ylabel("Error")
title("Problem 4A")
savefig("problem4A_better.png")

I_b_exact = 4.4889560612699568830
function I_b(x)
	return exp(cos((cos(x+1))^2+2*sin(4*x+1)))
end

y_0 = [error(newton_zero(I_b,0,pi,n), I_b_exact) for n in x]
y_1 = [error(newton_one(I_b,0,pi,n), I_b_exact) for n in x]
y_2 = [error(newton_two(I_b,0,pi,n), I_b_exact) for n in x]
using PyPlot
figure()
p = loglog(x,y_0,label="L0-th order")
loglog(x,y_1,label= L"1-st order")
loglog(x,y_2, label =L"2-nd order")
legend()
xlabel("N")
ylabel("Error")
title("Problem 4B")
savefig("problem4B_better.png")

I_c_exact = 6.6388149923287733132
function I_c(x)
	return tanh(x)/(.000000000000001+sqrt(abs(x-pi)))
end

#u-substitution
function I_c_left(u)
	return 2*tanh(pi-u^2)
end
function I_c_right(u)
	return 2*tanh(u^2+pi)
end

y_0 = [error(newton_zero(I_c_left,0,sqrt(pi),n/2.)+newton_zero(I_c_right,0,sqrt(pi),n/2.), I_c_exact) for n in x]
y_1 = [error(newton_one(I_c_left,0,sqrt(pi),n/2.)+newton_one(I_c_right,0,sqrt(pi),n/2.), I_c_exact) for n in x]
y_2 = [error(newton_two(I_c_left,0,sqrt(pi),n/2.)+newton_two(I_c_right,0,sqrt(pi),n/2.), I_c_exact) for n in x]
test_zero = [newton_zero(I_c_left,0,sqrt(pi),n/2.)+newton_zero(I_c_right,0,sqrt(pi),n/2.) for n in x]
using PyPlot
figure()
p = loglog(x,y_0,label=L"0-th order")
loglog(x,y_1,label= L"1-st order")
loglog(x,y_2, label =L"2-nd order")
legend()
xlabel("N")
ylabel("Error")
title("Problem 4C")
savefig("problem4C_better.png")

I_d_exact = 1.7981374998645790990
function I_d(x)
	return pi/(2*x)
end

y_0 = [error(newton_zero(I_d,1,pi,n), I_d_exact) for n in x]
y_1 = [error(newton_one(I_d,1, pi,n), I_d_exact) for n in x]
y_2 = [error(newton_two(I_d,1,pi,n), I_d_exact) for n in x]
using PyPlot
figure()
p = loglog(x,y_0,label=L"0-th order")
loglog(x,y_1,label= L"1-st order")
loglog(x,y_2, label =L"2-nd order")
legend()
xlabel("N")
ylabel("Error")
title("Problem 4D")
savefig("problem4D_better.png")

println("done")

