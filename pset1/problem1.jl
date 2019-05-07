println("hello")
function partial_sum(f, N)
	sum = 0
	for i in 1:N
		sum+=f(i)
	end
	return sum
end

function f(input)
	return 1/(input^4)
end
function error(psum)
	exact_soln = pi^4/90.
	num = abs(psum-exact_soln)
	return num/exact_soln
end
result = partial_sum(f, 1000)
println(result)
println(error(result))

x = [0,1,10,100,1000,10000,100000,1000000,10000000]
y = [error(partial_sum(f, i)) for i in x]
using PyPlot
#PyPlot.svg(true)


figure()
p = loglog(x,y, "o")
xlabel("N")
ylabel("Error")
title("Problem 1")
savefig("problem1.png")

