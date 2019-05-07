function kahan_sum(summand,N)
	sum=0
	c_last=0
	for i in 1:N
		y_i = summand(i, N)-c_last
		t_i = sum+y_i
		c_last = (t_i-sum)-y_i
		sum=t_i
	end
	return sum
end

function P_N(N)
	summand = pi/N
	S = 0.0
	for i=1:N
		S += summand
	end
	return S
end

function sum_direct(x, N)   
    S = 0.0
    for i=1:N
        S += x
    end
    return S
end

function sum_pairwise(x, N; base_threshold=128)
    if N < base_threshold
        return sum_direct(x, N)
    else
        N_half = Int64(floor(N/2))
        return sum_pairwise(x, N_half) + sum_pairwise(x, N - N_half)
    end
end

function summand(i,N)
	return pi/N
end

function relative_error(numerical, true_sol)
	numerator = abs(numerical-true_sol)
	return numerator/(1.0*true_sol)
end

true_sol = pi
N_vals = [10,100,1000,10000,100000,1000000,10000000,100000000]
kahan_err = [relative_error(kahan_sum(summand, n), true_sol) for n in N_vals]
pariwise_err = [relative_error(sum_pairwise(summand(0,n), n), true_sol) for n in N_vals]
regular_sum_err =  [relative_error(P_N(n), true_sol) for n in N_vals]

using PyPlot
figure()
loglog(N_vals,kahan_err,label=L"Kahan Sum")
loglog(N_vals,pariwise_err,label= L"Pairwise Sum")
loglog(N_vals, regular_sum_err, label =L"Naive Sum")
legend()
xlabel("Iterations N")
ylabel("Relative Error")
title("Problem 5 - Summations ")
savefig("problem5.png")
print("done")



