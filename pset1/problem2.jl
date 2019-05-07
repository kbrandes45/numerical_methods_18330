println("problem 2")

function P_N(N)
	summand = pi/N
	S = 0.0
	for i=1:N
		S += summand
	end
	return S
end

function error(psum)
	exact_soln = pi
	num = abs(psum-exact_soln)
	return num/exact_soln
end

x = [0,1,10,100,1000,10000,100000,1000000,10000000]
y = [error(P_N(i)) for i in x]

using PyPlot
figure()
p = loglog(x,y, "o")
xlabel("N")
ylabel("Error")
title("Problem 2")
savefig("problem2.png")
