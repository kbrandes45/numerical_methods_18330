using LinearAlgebra

function simpleMC(f, a, b, N)
	#returns <f>=sum/N instead of integral = V<f>
	D = length(a)
	S = 0.0
	for i=1:N
		x = rand(D)
		x = a .+ (b .- a).*x
		S+= f(x)
	end
	mean = S .* (1/N)
	std_dev = (S/N)^(.5)
	return std_dev #mean, std_dev
end


function coord_zone(N,p, a, b, f)
	D = length(a)
	samples = []
	for i=1:floor(Int, p*N)
		x = rand(D)
		x = a .+ (b .- a).*x
		push!(samples, x)
	end
	
	std_a = zeros(D) #i from 1,N
	count_a = zeros(D)
	std_b = zeros(D)
	count_b = zeros(D)
	for j=1:D
		#find all samples to match this direction
		a_per_zones = []
		b_per_zones = []
		for i=1:floor(Int, p*N)
			#iterate over each sample
			s = samples[i]
			if (check_in_A_zone(s, a, b, j))
				push!(a_per_zones, s)
			end
			if (check_in_B_zone(s, a, b, j))
				push!(b_per_zones, s)
			end
		end

		#now have list of all samples in the subvolume
		#compute stddev
		count_a[j] = length(a_per_zones)
		count_b[j] = length(b_per_zones)

		a_s = std_dev(a_per_zones, f)
		b_s = std_dev(b_per_zones, f)
		std_a[j] = a_s
		std_b[j] = b_s
	end

	#find minimal
	index=0
	best_so_far = 100000000
	for i in 1:D
		t = std_a[i]+std_b[i]
		if (t<best_so_far)
			index = i
			best_so_far = t
		end
	end

	new_a1 = []
	new_b1 = []
	new_a2 = []
	new_b2 = []
	for i=1:D
		if (i==index)
			ai=a[index]
			bi=b[index]
			midpt = (ai+bi)/2.0
			push!(new_a1, ai)
			push!(new_b1, midpt)

			push!(new_a2, midpt)
			push!(new_b2, bi)
		else
			push!(new_a1, a[i])
			push!(new_b1, b[i])
			push!(new_a2, a[i])
			push!(new_b2, b[i])
		end
	end

	N_A1 = count_a[index]
	N_B2 = count_b[index]
	std_A1 = std_a[index]
	std_B2 = std_b[index]

	if ((std_A1 == 0)||(std_B2==0))
		N_A1_final = .5*(1-p)*N+N_A1
		N_B2_final = .5*(1-p)*N+N_B2
	else
		NA_new = ((std_A1)/(std_A1+std_B2))*(1-p)*N
		N_A1_final = count_a[index]+NA_new
		N_B2_final = count_b[index] + ((1-p)*N-NA_new)
	end
	return (N_A1_final, new_a1, new_b1), (N_B2_final, new_a2, new_b2)
end


function recurse_strat_samp(f, a, b, N, p)
	if (N<=100)
		current_std = simpleMC(f, a, b, N) 
		return current_std #current_mean, current_std
	else

		(N_A1_final, new_a1, new_b1), (N_B2_final, new_a2, new_b2) = coord_zone(N,p, a, b, f)
		#tot_mean = .5*(current_meanA+current_meanB)
		#tot_std_dev = .5*(current_stdA^2+current_stdB^2)^(.5)
		return .5*(recurse_strat_samp(f, new_a1, new_b1, N_A1_final, p)^2+recurse_strat_samp(f, new_a2, new_b2, N_B2_final, p)^2)
	end
end

function std_dev(samples, f)
	#i think N is same as length(samples)
	if (length(samples) <= 0)
		return 0.0
	end
	su = 0.0
	for s=1:length(samples)
		samp = samples[s]
		su+=f(samp)
	end
	return (su/length(samples))^(.5)
end


function check_in_A_zone(x, a, b, i)
	x_si = x[i]
	if (a[i]<= x_si)
		if (x_si< (a[i]+b[i])/2.0)
			return true
		end
	end
	return false
end

function check_in_B_zone(x, a, b, i)
	x_si = x[i]
	if (((a[i]+b[i])/2.0)<= x_si)
		if (x_si<= b[i])
			return true
		end
	end
	return false
end



function f(x)
    if (norm(x)<1.0)
        return 1.0
    else
        return 0.0
    end
end

a= [0.0,0.0]
b=[1.0,1.0]
N=10000
p = .01
std_devs= recurse_strat_samp(f, a, b, N, p)
println("DONE")
println(std_devs)


Ns = [1000, 10000, 100000, 1000000]
std_per_n = []
for n in Ns
	current_std_devs = []
	for i=1:100
		println("Iteration ",i)
		std_devs= recurse_strat_samp(f, a, b, n, p)
		push!(current_std_devs, std_devs)
	end
	push!(std_per_n, current_std_devs)
end


using PyPlot
figure()
#exact=.5
#axhline(exact, color="k", ls="--", label="exact value")
semilogy(std_per_n[1], "o", label=L"N=1000")
semilogy(std_per_n[2], "o", label=L"N=10000")
semilogy(std_per_n[3], "o", label=L"N=100000")
println(std_per_n[3])
semilogy(std_per_n[4], "o", label=L"N=1000000")

ylabel("Monte Carlo Approximation Std Dev")
xlabel("Sample index")

legend()

title("Problem 3")
savefig("problem3_logscale.png")
print("done")

Ns = [1000, 10000, 100000, 1000000]
std_per_n_simple = []
for n in Ns
	current_std_devs = []
	for i=1:100
		println("Iteration ",i)
		std_devs= simpleMC(f, a, b, N)
		push!(current_std_devs, std_devs)
	end
	push!(std_per_n_simple, current_std_devs)
end


function get_average(s)
	return sum(s)/length(s)
end

avg = [get_average(s) for s in std_per_n]
avgsimple = [get_average(s) for s in std_per_n_simple]
figure()
loglog(Ns, avg, "o-", label="Recursive Stratified MC")
loglog(Ns, avgsimple, "o-", label="Simple MC")

legend()
ylabel("Monte Carlo Approximation Std Dev (Averaged)")
xlabel("N")


title("Problem 3 Averaged Std Devs")
savefig("problem3_avgstd.png")
print("done")