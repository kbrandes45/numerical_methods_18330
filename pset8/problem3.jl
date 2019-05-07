using FFTW
using LinearAlgebra 

function solver(c,b)
	#fft(inv(diag(fft(c))))*ifft(b)
	b_updated = ifft(b)
	evals = fft(c)
	inv_evals = 1 ./ evals
	diag_evals = Diagonal(inv_evals)
	bbar = diag_evals*b_updated
	soln = fft(bbar)
	return soln
end

#Test problem for n=4
n=4
b = rand(Int, n)
c = rand(Int, n)

C = zeros(n, n)
for i = 1:n
   for j = 1:n-i+1
       C[i,i+j-1] = c[j]
   end
   for j = n-i+2:n
       C[i, j-(n-i+1)] = c[j]
   end
end
backslash = C\b
println(backslash)
s=solver(c,b)
println(s)


ns = [10, 100, 1000, 10000]
ns = [10, 20, 50, 100, 250, 500, 1000, 2000, 5000, 10000]

back_time = []
solve_time = []
for i in 1:length(ns)
	n = ns[i]

	b = rand(Int, n)
	c = rand(Int, n)

	C = zeros(n, n)
	for i = 1:n
	   for j = 1:n-i+1
	       C[i,i+j-1] = c[j]
	   end
	   for j = n-i+2:n
	       C[i, j-(n-i+1)] = c[j]
	   end
	end

	backslash_time = @timed C\b
	b_t = backslash_time[2]
	#println(backslash_time)
	println(b_t)
	push!(back_time, b_t)
	solver_time = @timed solver(c,b)
	s_t = solver_time[2]
	push!(solve_time, s_t)
end

using PyPlot
figure()
print(back_time)
semilogx(ns, back_time , "o-", label="BackSlash Operator")
semilogx(ns, solve_time, "o-", label="FFT Circulant System Solver")
xlabel("N")
ylabel("Time")

legend()

title("Problem 3")
savefig("problem3_times.png")
print("done")


#solving huge system
n=10000000
b = rand(Int, n)
c = rand(Int, n)

#C = zeros(n, n)
#for i = 1:n
#   for j = 1:n-i+1
#       C[i,i+j-1] = c[j]
#   end
#   for j = n-i+2:n
#       C[i, j-(n-i+1)] = c[j]
#   end
#end
#print(C\b)
println(solver(c,b))

