print("starting")

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

function adaptive(f,a,b,quadrature,N,epsilon)
	q_n = quadrature(f,a,b,N)
	q_two_n = quadrature(f,a,b,2*N)
	error_est = abs(q_two_n-q_n)
	if (error_est<epsilon) 
		return q_two_n
	else
		return adaptive(f,a,(a+b)/2.0,quadrature,N,epsilon)+adaptive(f,(a+b)/2.0,b,quadrature,N,epsilon)
	end
end

function test_function(x)
	return (x^10)*exp(4*x^3-3*x^4)
end

function relative_error(numerical, true_sol)
	numerator = abs(numerical-true_sol)
	return numerator/(1.0*true_sol)
end


#Define interval and N
N=1
a=0
b=3
true_soln = newton_two(test_function, a, b, 10^5)
print("here")
epsilon = [.1,.01,.001,.0001,.00001,.000001, .0000001]

print(adaptive(test_function,a,b,newton_zero,N,.1))

y_0 = [relative_error(adaptive(test_function,a,b,newton_zero,N,epi), true_soln) for epi in epsilon]
print("done with y_0")
y_1 = [relative_error(adaptive(test_function,a,b,newton_one,N,epi), true_soln) for epi in epsilon]
print("done with y_1")
y_2 = [relative_error(adaptive(test_function,a,b,newton_two,N,epi), true_soln) for epi in epsilon]

using PyPlot
figure()
loglog(epsilon,y_0,label=L"0-th order")
loglog(epsilon,y_1,label= L"1-st order")
loglog(epsilon,y_2, label = L"2-nd order")
legend()
xlabel("Epsilon")
ylabel("Relative Error")
title("Problem 1 - Adaptive Quadrature")
savefig("problem1.png")
print("done")