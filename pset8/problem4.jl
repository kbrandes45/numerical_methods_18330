using FFTW
using LinearAlgebra

#To compare, do n=3
n=3
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
println(eigvals(C))
println(fft(c))


ns = [10, 20, 50, 100, 250, 500, 1000]
f_time = []
e_time = []
for i in 1:length(ns)
	n = ns[i]

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

	fft_time = @timed fft(c)
	res = fft_time[1]
	f_t = fft_time[2]
	push!(f_time, f_t)

	eval_time = @timed eigvals(C)
	res1 = eval_time[1]
	e_t = eval_time[2]
	push!(e_time, e_t)
end

using PyPlot
figure()
semilogx(ns, e_time , "o-", label="LinearAlgebra eigvals Function")
semilogx(ns, f_time, "o-", label="FFT Eigenvalue Solver")
xlabel("N")
ylabel("Time")

legend()

title("Problem 3")
savefig("problem4_times.png")
print("done")
