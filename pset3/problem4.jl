function split_N(N)
	if N%2==0
		return [N/2, N/2]
	else
		return [N/2-.5, N/2+.5]
	end
end

function simpsons(f, a, b, N)
	delta = (b-a)/(1.0*N)
	h = delta/2.0
	sum = h*f(a)/3+h*f(b)/3
	for n in 1:N-1
		sum+=(h/3)*(2*f(a+2*n*h)+4*f(a+(2*n-1)*h))
	end
	sum+= 4*h*f(a+(2*N-1)*h)/3
	return sum
end

function recursive_simpsons(f,a,b,N)
	if N<100
		return simpsons(f,a,b,N)
	else
		Ns = split_N(N)
		midpt = (a+b)/2.0
		return recursive_simpsons(f,a,midpt,Ns[1])+recursive_simpsons(f,midpt,b,Ns[2])
	end
end

function f(x)
	return (Base.log(x))^3 #natural log
end

function relative_error(numerical, true_sol)
	numerator = abs(numerical-true_sol)
	return numerator/(1.0*true_sol)
end

a=1
b=2
true_solution = 0.10109738718799412444
x = [1,10,100,1000,10000,100000,1000000,10000000,100000000]
y_simpson_only = [simpsons(f,a,b,i) for i in x]
y_recursive_simpsons= [recursive_simpsons(f,a,b,i) for i in x]
println(y_simpson_only)
println(y_recursive_simpsons)
y_simp_err = [relative_error(i, true_solution) for i in y_simpson_only]
y_rec_err = [relative_error(i, true_solution) for i in y_recursive_simpsons]

using PyPlot
figure()
loglog(x,y_simp_err,label=L"Regular Simpsons Method")
loglog(x,y_rec_err,label= L"Recursive Simpsons Method")
legend()
xlabel("Iterations N")
ylabel("Relative Error")
title("Problem 4 - Recursive Simpsons ")
savefig("problem4.png")
print("done")





