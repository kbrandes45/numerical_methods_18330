# Jenkins Traub RPOLY
# Representation for a polynomial p(x) = a_0+a_1x+a_2x^2+...+a_{n}x^{n} will be a array of the coeeficients
# So p(x) = [a_n, a_{n-1}, ... , a_2, a_1, a_0]
#FIRST COEFF = LEADING COEFF = HIGHEST ORDER TERM = INDEX ONE OF P_COEFFS
# Note that all the coefficients should be in the list - if they do have a value, they should be set to 0.

function rpoly_stage1(p_coeffs, num_iters, epsilon)
	#make zero-shift polynomial sequences; approximately 5 iterations

	#Take the derivative
	powers = [length(p_coeffs)-i for i in 1:length(p_coeffs)]
	k_0 =  powers .* p_coeffs
	k_0 = k_0[1:length(k_0)-1]

	#normalize the derivative function
	K_polys = k_0 ./ k_0[1]
	counter = 0

	#Iterate through fixed number of 0-shifts; terminate if root is found
	for i in 1:num_iters
		K_polys, x, root_found = next_step(p_coeffs, K_polys, 0, epsilon, false)
		counter+=1
		if (root_found)
			break
		end
	end
	return K_polys, counter
end

function rpoly_stage2(p_coeffs, K_curr, lin_const, num_iters, epsilon)
	#fixed shift stage 2 

	t_next = t = Inf
	counter = 0
	
	#normalize coefficients and operate off of this
	if (K_curr[1] != 1)
		K_polys = K_curr ./ K_curr[1] 
	else
		K_polys = K_curr
	end

	#Iteratively compute fixed-shifts until a sufficient root approximation is found (or maximum number of iterations exceeded)
	root_found = false
	while (true)
		t, t_prev = t_next, t
		K_polys, t_next, root_found = next_step(p_coeffs, K_polys, lin_const, epsilon, true)
		if (root_found)
			break
		end

		#Termination conditions combaring the difference in root approximations
		counter += 1
		cond1 = abs(t-t_prev) <= .5*abs(t_prev) #is this a strictly less than or a less than
		cond2 = abs(t_next-t) <= .5*abs(t)

		if ((cond1 && cond2) || (counter > num_iters))
			break
		end
	end

	return K_polys, counter, root_found

end


function rpoly_stage3(p_coeffs, K_curr, lin_const, num_iters, epsilon)
	#stage 3 variable shifts

	#normalize coefficients and operate off of this
	if (K_curr[1] != 1)
		K_curr = K_curr ./ K_curr[1] 
	end

	#Set initial estimate s_L for the root
	sl = lin_const - synth_div_eval(p_coeffs, lin_const)[2]/synth_div_eval(K_curr, lin_const)[2]
	sl_prev = Inf
	counter = 0 

	K_polys = K_curr

	#Iterate until the absolute difference root approximations for 2 iterations is < epsilon
	#Or also terminate if exceed maximum number of allowed iterations
	while ((abs(sl - sl_prev) > epsilon) && (counter < num_iters))
		p_prime, p_eval_at_sl = synth_div_eval(p_coeffs, sl)
		k_next, k_curr_at_sl = synth_div_eval(K_polys, sl)
		k_next = pushfirst!(k_next, 0)

		#check if we find a root by seeing if evaluated polynomial is 0
		if (abs(p_eval_at_sl) < epsilon)
			return K_polys[length(K_polys)], sl, counter, true
		end	
		
		k_next_iter = p_prime .- (p_eval_at_sl / k_curr_at_sl) .* k_next
		counter+=1
		sl, sl_prev = (sl - p_eval_at_sl / synth_div_eval(k_next_iter, sl)[2]), sl
		K_polys = k_next_iter
	end
	if (counter > num_iters)
		#Failed to find a root for this iteration, returning false
		return 0,0,0,false
	end
	return K_polys[length(K_polys)], sl, counter, true
end

function newton(start, func_coeffs, deriv_coeffs, epsilon, num_iters)
	#Iterative newton raphson algorithm that terminates when aboslute error difference is < epsilon
	x_curr = start
	x_next = x_curr - evaluate_poly(func_coeffs, x_curr)/evaluate_poly(deriv_coeffs, x_curr)
	counter = 0
	while ((abs(x_next - x_curr)>epsilon) && (counter < num_iters))
		x_curr = x_next
		x_next = x_curr - evaluate_poly(func_coeffs, x_curr)/evaluate_poly(deriv_coeffs, x_curr)
		counter += 1
	end
	return x_curr
end

function jenkins_traub_one_iter(p_coeffs, num_iterations, epsilon)
	#Compute the smallest root for the current coefficients list using the 3 stage Jenkins Traub alg

	#Check if the constant coeff is a zero --> implies the root is zero
	if ((length(p_coeffs) > 0) && (p_coeffs[length(p_coeffs)]==0))
		return p_coeffs, 0, 0 
	end

	#check if the current poly is only a constant.
	if (length(p_coeffs) < 2)
		return 
	end

	#normalize coefficients and operate off of this
	if (p_coeffs[1] != 1)
		p_coeffs = p_coeffs ./ p_coeffs[1] 
	end

	if (length(p_coeffs)==2)
		#linear function,so root is just the constant. return none for coefs since nothing is left
		return 0.0, -p_coeffs[length(p_coeffs)], 0 
	end
	total_func_evals = 0
	k_curr, num_iters = rpoly_stage1(p_coeffs, 5, epsilon)
	total_func_evals += 2*num_iters #next_step has 2 func evals
	println("Stage 1 with ",num_iters)

	#get a guess for initial linear constant s; approximate using standard method of a cauchy polynomial + netwons iter
	# new polynomial has coeffs that are the moduli of p_coeffs, but the last root is -1.
	all_pos = [abs(i) for i in p_coeffs] #TODO: need an element wise absolute value
	all_pos[length(all_pos)] = -all_pos[length(all_pos)]

	#Take the derivative
	powers = [length(all_pos)-i for i in 1:length(all_pos)]
	deriv =  powers .* all_pos
	deriv = deriv[1:length(deriv)-1]

	beta = newton(1.0, all_pos, deriv, .01, 500)

	while (true)
		phi_rand = 2.0*pi*rand(1)[1] 
		sl = cos(phi_rand)*beta

		k_curr, num_iters2, root_found = rpoly_stage2(p_coeffs, k_curr, sl, num_iterations, epsilon)
		println("Done with stage 2 at ",num_iters2," iterations")
		total_func_evals += num_iters2*2 #next_step has 2 evals

		k_curr, sl, num_iters3, success =  rpoly_stage3(p_coeffs, k_curr, sl, num_iterations, epsilon)
		total_func_evals += num_iters3*3 #next_step and one extra eval
		println("Stage 3 finsihed with total iterations ",num_iters3)
		if (!success)
			#didn't converge, try increasing maxiterations 
			num_iterations = 2*num_iterations
		else
			return k_curr, sl, total_func_evals
		end
	end

end

function get_all_roots_JT(p_coeffs, epsilon, num_iterations)
	#Return the total set of coefficients for a polynomial and the root

	#Remove any leading zero coefficients since they are simply extra terms and will mess up normalization
	while ((length(p_coeffs) > 0) && (p_coeffs[1]==0))
		p_coeffs = p_coeffs[2:length(p_coeffs)]
	end

	if (length(p_coeffs) < 2)
		# function was just a constant, no roots
		return 
	end

	roots = []
	total_function_evals = 0
	#Iteratively run Jenkins-Traub algorithm
	for i in 1:length(p_coeffs)-1
		upd_coeffs, s, num_f_evals = jenkins_traub_one_iter(p_coeffs, num_iterations, epsilon)

		#update the polynomial by dividing out the linear term that factors in the root
		p_coeffs, eval = synth_div_eval(p_coeffs, s)

		#increase count on number of times the polynomial func is evaluated
		total_function_evals += num_f_evals+1

		#store the root found for the iteration
		roots = push!(roots, s)
	end
	return roots, total_function_evals
end

function evaluate_poly(p_coeffs, val)
	#Stand alone polynomial evaluation method (no synthetic division)
	highest_power = length(p_coeffs)
	sum = 0
	for i in 1:highest_power
		power = highest_power-i
		sum+= p_coeffs[i]*val^power
	end
	return sum
end

function synth_div_eval(p_coeffs, lin_const)
	#Assumes leading coeff is 1.0 Decrease length of list by 1 (since dividing by x+/-c)


	#To only divide the entire function by x, then set lin_const = 0. This involves only truncating the function. Evaluation is merely the constant term
	if (lin_const == 0)
		return p_coeffs[1:length(p_coeffs)-1], p_coeffs[length(p_coeffs)]
	end

	#In all other cases, perform regular synthetic division
	syn_div = zeros(length(p_coeffs)-1)
	syn_div[1] = p_coeffs[1]
	for i=2:length(p_coeffs)-1
		syn_div[i] = p_coeffs[i] + syn_div[i-1]*lin_const
	end

	#Evaluate linear portion of polynomial -- taking advantage of remainder of the poly division
	eval = p_coeffs[length(p_coeffs)] + syn_div[length(syn_div)] .* lin_const

	return syn_div, eval
end

function next_step(p_coeffs, K_current, const_to_eval, epsilon, generate_t)
	#Update polynomial then evaluate at const; epsilon = threshold for roots

	p_prime, p_eval_at_c = synth_div_eval(p_coeffs, const_to_eval)
	k_next, k_curr_at_c = synth_div_eval(K_current, const_to_eval)

	root_found = false

	#check if we find a root by seeing if evaluated polynomial is 0
	if (abs(p_eval_at_c) < epsilon)
		root_found = True
	end

	# If evaluated shifted polynomial is too small, then increase it slightly so not dividing by 0
	if (abs(k_curr_at_c) < epsilon)
		k_curr_at_c += epsilon/100.0
	end

	t = 0.0
	if (generate_t == true)
		#update root approximation for the new shifted polynomial (not needed for all stages)
		t = const_to_eval - p_eval_at_c/k_curr_at_c
	end

	#shift by zero once to represent the divide by 0
	k_deriv = pushfirst!(k_next, 0.0)

	#update the polynomial
	final = p_prime .- (p_eval_at_c / k_curr_at_c) .* k_deriv

	return final, t, root_found
end


############# Total Testing ###################
#polynomial = x^2-5x+6 = (x-2)(x-3) so roots = 2,3
roots, fc = get_all_roots_JT([1.0,-5.0,6.0], .00001, 100)
println(roots)


#polynomial= (x+1) so roots = -1
roots, fc = get_all_roots_JT([1.0,1.0], .0000001, 100)
println(roots)

#polynomial = (2x-3) so roots = 3/2
roots, fc = get_all_roots_JT([2.0, -3.0], .0000001, 100)
println(roots)

#polynomial = (x^2-4x) so roots = 0, 4
roots, fc = get_all_roots_JT([1.0, -4.0, 0.0], .0000001, 100)
println(roots)

#polynomial = (x-1)(x-2)(x+3)(x+4)(x-5)(x+6) --> roots = 1,2,-3,-4,5,-6
# [1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0]
roots, fc = get_all_roots_JT([1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0], .00000001, 1000)
println(roots)


#polynomial = (x-1)(x-2)(x+3)(x+4)(x-5)(x+6)(x+45)(x-2.2)(x+3.77)
roots,fc = get_all_roots_JT([1.0, 0.0, 51.57, 262.206, -1747.26, 10650.8, 13582.4, 99964.6, -70003.7, -300186.0, 268726.0], .00000001, 100000)
println(roots)


###########To evaluate versus newton rapheson
function rel_err(approx, exact)
	return abs(abs(approx-exact)/exact)
end

function newton_raphson(f, f_prime, x0, N)
    """ Implement Newton's method for finding a
    zero of f(x) from the starting guess x0.
    """
    x = x0
    for i=1:N
        x = x .- f(x)./f_prime(x)
    end
    
    return x
end

function test_func(x)
	return evaluate_poly([1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0], x)
end

function der_test(x)
	p_coeffs  = [1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0]
	powers = [length(p_coeffs)-i for i in 1:length(p_coeffs)]
	k_0 =  powers .* p_coeffs
	k_0 = k_0[1:length(k_0)-1]
	return evaluate_poly(k_0, x)
end

#Newton Raphson will find smallest root, so will inner JT
Ns = [1,2,3,4,5,6,7,8,9,10]
guess = 10.0  
n_r_roots = [newton_raphson(test_func, der_test, guess, n) for n in Ns]
println(n_r_roots)
p_coeffs = [1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0]
j_t_roots = [jenkins_traub_one_iter(p_coeffs, n, 1e-10)[2] for n in Ns]
println(j_t_roots)

n_r_err = [rel_err(n_r_roots[i], 5.0) for i in 1:length(Ns)]
j_t_err = [rel_err(j_t_roots[i], 1.0) for i in 1:length(Ns)]
using PyPlot
figure()
plot(Ns, n_r_err , "o-", label="Newton Raphson Method")
plot(Ns, j_t_err, "o-", label="Jenkins Traub for 1 (Real) Root")
xlabel("N")
ylabel("Relative Error")

legend()

title("Error Convergence Rates (Small Polynomial)")
savefig("jenkins_error_rates.png")
print("done")

#####TEST for a medium polynomial (degree = 15)
function test_func(x)
	return evaluate_poly([280.8, 3467.81, 781.739, -128636.0, -347272.0, 1.1595e6, 4.64429e6, -2.54983e6, -2.20585e7, -9.42948e6, 3.89638e7, 4.0433e7, -9.45315e6, -2.83132e7, -1.175e7, -1.17e6], x)
end

function der_test(x)
	p_coeffs  = [280.8, 3467.81, 781.739, -128636.0, -347272.0, 1.1595e6, 4.64429e6, -2.54983e6, -2.20585e7, -9.42948e6, 3.89638e7, 4.0433e7, -9.45315e6, -2.83132e7, -1.175e7, -1.17e6]
	powers = [length(p_coeffs)-i for i in 1:length(p_coeffs)]
	k_0 =  powers .* p_coeffs
	k_0 = k_0[1:length(k_0)-1]
	return evaluate_poly(k_0, x)
end


big_boi = [280.8, 3467.81, 781.739, -128636.0, -347272.0, 1.1595e6, 4.64429e6, -2.54983e6, -2.20585e7, -9.42948e6, 3.89638e7, 4.0433e7, -9.45315e6, -2.83132e7, -1.175e7, -1.17e6]
Ns = [1,2,3,4,5,6,7,8,9,10,11,12,13]
guess = 0#-7
n_r_roots = [newton_raphson(test_func, der_test, guess, n) for n in Ns]
j_t_roots = [jenkins_traub_one_iter(big_boi, n, 1e-10)[2] for n in Ns]
println(n_r_roots)
println(j_t_roots)
println(get_all_roots_JT(big_boi, 1e-10, 10))
exact = -.1489

n_r_err = [rel_err(n_r_roots[i], exact) for i in 1:length(Ns)]
j_t_err = [rel_err(j_t_roots[i], exact) for i in 1:length(Ns)]

using PyPlot
figure()
semilogy(Ns, n_r_err , "o-", label="Newton Raphson Method")
semilogy(Ns, j_t_err, "o-", label="Jenkins Traub for 1 (Real) Root")
xlabel("N")
ylabel("Relative Error")

legend()

title("Error Convergence Rates (Medium Polynomial)")
savefig("jenkins_error_rates_numone.png")


###### TEST for a large polynomial (big spread of roots) = 26th degree
#expand (x-10)(x+15)(2x-13.99)(4x+25)(x+56)(-3x-1)(-30x+5)(21x+16)(x-90)(x-1)(x-2)(x+3)(x+4)(x-5)(x+6)(2x+3)( x-2.3)(-3x-3)(2x+.3)(-x-.66)(x+35)(2x+1.2)(x-7.2)(x-72)(x+28)(x- .1)
biggg_coeffs = [362880.0, -1.39732e7, -3.17741e9, 2.52633e10, 8.36824e12, 1.54068e14, -1.52172e15, -3.86931e16, 9.95976e16, 3.20459e18, -1.33643e18, 1.15358e20, -8.77691e19, 1.7984e21, 2.75588e21, -1.00243e22, -2.07459e22, 1.38201e22, 5.2646e22, 2.23018e22, -2.74017e22, -2.8284e22, -7.53574e21, 4.93807e20, 3.74285e20, 6.29872e18, -3.80548e18]
function test_func(x)
	return evaluate_poly([362880.0, -1.39732e7, -3.17741e9, 2.52633e10, 8.36824e12, 1.54068e14, -1.52172e15, -3.86931e16, 9.95976e16, 3.20459e18, -1.33643e18, 1.15358e20, -8.77691e19, 1.7984e21, 2.75588e21, -1.00243e22, -2.07459e22, 1.38201e22, 5.2646e22, 2.23018e22, -2.74017e22, -2.8284e22, -7.53574e21, 4.93807e20, 3.74285e20, 6.29872e18, -3.80548e18], x)
end

function der_test(x)
	p_coeffs  = [362880.0, -1.39732e7, -3.17741e9, 2.52633e10, 8.36824e12, 1.54068e14, -1.52172e15, -3.86931e16, 9.95976e16, 3.20459e18, -1.33643e18, 1.15358e20, -8.77691e19, 1.7984e21, 2.75588e21, -1.00243e22, -2.07459e22, 1.38201e22, 5.2646e22, 2.23018e22, -2.74017e22, -2.8284e22, -7.53574e21, 4.93807e20, 3.74285e20, 6.29872e18, -3.80548e18]
	powers = [length(p_coeffs)-i for i in 1:length(p_coeffs)]
	k_0 =  powers .* p_coeffs
	k_0 = k_0[1:length(k_0)-1]
	return evaluate_poly(k_0, x)
end

Ns = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
guess = 0
n_r_roots = [newton_raphson(test_func, der_test, guess, n) for n in Ns]
j_t_roots = [jenkins_traub_one_iter(biggg_coeffs, n, 1e-10)[2] for n in Ns]

#Double exact values because Jenkins Traub will converge to either of these depending on the random initial shift values in stage 2
exact = .1
exact2 = -.15

function special_rel(approx, exact1, exact2)
	if (abs(approx-exact1) < abs(approx-exact2))
		return rel_err(approx, exact1)
	else
		return rel_err(approx, exact2)
	end
end

n_r_err = [rel_err(n_r_roots[i], .1666666666666) for i in 1:length(Ns)]
j_t_err = [special_rel(j_t_roots[i], exact, exact2) for i in 1:length(Ns)]

using PyPlot
figure()
semilogy(Ns, n_r_err , "o-", label="Newton Raphson Method")
semilogy(Ns, j_t_err, "o-", label="Jenkins Traub for 1 (Real) Root")
xlabel("N")
ylabel("Relative Error")

legend()

title("Error Convergence Rates (Large Polynomial)")
savefig("jenkins_error_rates_numtwo.png")

function newton_raphson(fc, fc_prime, x0, N)
    x = x0
    for i=1:N
        x = x .- evaluate_poly(fc, x)./evaluate_poly(fc_prime, x)
    end
    return x
end

Ns = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37, 38,39,40,41,42,43,44,45,46,47,48]
guess = 10000
c = [1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0]
powers = [length(c)-i for i in 1:length(c)]
cprime =  powers .* c
cprime = cprime[1:length(cprime)-1]
n_r_roots = [newton_raphson(c, cprime, guess, n) for n in Ns]
n_r_roots_g = [newton_raphson(c, cprime, 10.0, n) for n in Ns]
exact = 5.0

n_r_err = [rel_err(n_r_roots[i], exact) for i in 1:length(Ns)]
n_r_err_g = [rel_err(n_r_roots_g[i], exact) for i in 1:length(Ns)]

using PyPlot
figure()
semilogy(Ns, n_r_err , "o-", label="Bad Initial Guess")
semilogy(Ns, n_r_err_g , "o-", label="Good Initial Guess")

xlabel("N")
ylabel("Relative Error")
legend()

title("Newton Raphson Convergence with Bad Initial Guess")
savefig("NR_badguess.png")
print("done")


function iterative_newton_raphson(c, N)
	#Newton Raphson for polynomials of all degrees
	powers = [length(c)-i for i in 1:length(c)]
	cprime =  powers .* c
	cprime = cprime[1:length(cprime)-1]
	roots = []
	for i in 1:length(c)-1
		r = newton_raphson(c, cprime, 0.0, N)
		push!(roots, r)
		c = synth_div_eval(c, r)[1]
		powers = [length(c)-i for i in 1:length(c)]
		cprime =  powers .* c
		cprime = cprime[1:length(cprime)-1]
	end
	return roots
end

c = [1.0, 5.0, -33.0, -149.0, 212.0, 684.0, -720.0]

roots_NR = iterative_newton_raphson(c, 100)
roots_JT = get_all_roots_JT(c, 1e-10, 10)[1]
println(roots_NR)
println(roots_JT)



#Compare time to solve with number of roots
biggg_coeffs = [362880.0, -1.39732e7, -3.17741e9, 2.52633e10, 8.36824e12, 1.54068e14, -1.52172e15, -3.86931e16, 9.95976e16, 3.20459e18, -1.33643e18, 1.15358e20, -8.77691e19, 1.7984e21, 2.75588e21, -1.00243e22, -2.07459e22, 1.38201e22, 5.2646e22, 2.23018e22, -2.74017e22, -2.8284e22, -7.53574e21, 4.93807e20, 3.74285e20, 6.29872e18, -3.80548e18]
#25 roots

function take_deriv(c)
	powers = [length(c)-i for i in 1:length(c)]
	cprime =  powers .* c
	cprime = cprime[1:length(cprime)-1]
	return cprime
end
function generate_all_degree_polys(biggg_coeffs)
	#Generate polynomials with real valued coefficients by iteratively taking the derivative and storing that as a standalone polynomial 
	#Creates N-2 (where N=degree of initial coefficients list) total polynomials
	many_poly = [biggg_coeffs]
	for i in 1:length(many_poly[1])-2
		c = take_deriv(many_poly[length(many_poly)])
		many_poly = push!(many_poly, c)
	end
	return many_poly
end
many_poly = generate_all_degree_polys(biggg_coeffs)

solve_times = []
poly_length = []
function get_time(c, solve_times, poly_length)
	res = @timed get_all_roots_JT(c, 1e-10, 10)
	time = res[2]
	solve_times=push!(solve_times, time)
	poly_length=push!(poly_length, length(c))
end

function get_time2(c, solve_times, poly_length)
	res = @timed iterative_newton_raphson(c, 100)
	time = res[2]
	solve_times=push!(solve_times, time)
	poly_length=push!(poly_length, length(c))
end

for i in 1:length(many_poly)
	get_time(many_poly[length(many_poly)-i+1], solve_times, poly_length)
	println(length(many_poly[length(many_poly)-i+1]))
end

solve_times_NR = []
poly_length_NR = []
for i in 1:length(many_poly)
	get_time2(many_poly[length(many_poly)-i+1], solve_times_NR, poly_length_NR)
end
using PyPlot
figure()
plot(poly_length, solve_times, "o-", label = "Jenkins Traub Method")
plot(poly_length_NR, solve_times_NR , "o-", label="Newton Raphson Method")
xlabel("Polynomial Degree/Order")
ylabel("Solve Times [seconds]")
legend()

title("Solve Times for Computing All Roots on Varying Degree Polynomials")
savefig("solve_times.png")
print("done")


c = [280.8, 3467.81, 781.739, -128636.0, -347272.0, 1.1595e6, 4.64429e6, -2.54983e6, -2.20585e7, -9.42948e6, 3.89638e7, 4.0433e7, -9.45315e6, -2.83132e7, -1.175e7, -1.17e6]


#COUNTING FUNCTION EVALS
function newton_raphson(fc, fc_prime, x0, N)
    x = x0
    counter = 0
    x_last = 1000000
    while (abs(x-x_last)>1e-10)
    	x_last = x
    	counter+=2
        x = x .- evaluate_poly(fc, x)./evaluate_poly(fc_prime, x)
        if (counter > 5000)
        	println("exceeding counter max")
        	break
        end
    end
    return x, counter
end

function iterative_newton_raphson_count_evals(c, N)
	powers = [length(c)-i for i in 1:length(c)]
	cprime =  powers .* c
	cprime = cprime[1:length(cprime)-1]
	roots = []
	counter = 0
	for i in 1:length(c)-1
		r, count = newton_raphson(c, cprime, 0.0, N)
		push!(roots, r)
		counter+=count
		c = synth_div_eval(c, r)[1]
		powers = [length(c)-i for i in 1:length(c)]
		cprime =  powers .* c
		cprime = cprime[1:length(cprime)-1]
		counter += 1 #for the synthetic division
	end
	return roots, counter
end

many_poly = generate_all_degree_polys(biggg_coeffs)

function count_evals(many_poly)
	JT_evals = []
	NR_evals = []
	poly_length = []
	for i in 1:length(many_poly)
		jt_e = get_all_roots_JT(many_poly[length(many_poly)-i+1], 1e-10, 10)[2]
		JT_evals = push!(JT_evals, jt_e)
		nr_e = iterative_newton_raphson_count_evals(many_poly[length(many_poly)-i+1], 100)[2]
		NR_evals = push!(NR_evals, nr_e)
		poly_length = push!(poly_length, length(many_poly[length(many_poly)-i+1]))
	end
	return JT_evals, NR_evals, poly_length
end

JT_evals, NR_evals, poly_length = count_evals(many_poly)

using PyPlot
figure()
semilogy(poly_length, JT_evals, "o-", label = "Jenkins Traub Method")
semilogy(poly_length, NR_evals , "o-", label="Newton Raphson Method")
xlabel("Polynomial Degree/Order")
ylabel("Function Evaluations")
legend()

title("Number of Function Evaluations for Varying Degree Polynomials")
savefig("func_evals.png")
print("done")
