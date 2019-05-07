
function richardson_extrapolate(F, Δ, t, p)
    """ Take a numerical method F(Δ) of order p and return its Richardson
    extrapolated version.
    """
    over_t, over_t_calls = F(Δ/t)
    f_delt, f_delt_calls = F(Δ)
    return (t^p * over_t - f_delt)/(t^p - 1), over_t_calls+f_delt_calls
end

function newton_zero(f,a,b, delta)
	N = (b-a)/(1.0*delta)
	sum=0.0
	count =0.0
	for i in 1:N-1
		sum+= delta*f(a+i*delta)
		count+=1
	end
	return sum, count
end 

function newton_one(f,a,b, delta)
	N = (b-a)/(1.0*delta)
	sum=.5*delta*(f(a)+f(b))
	count = 2
	for i in 1:N-1
		sum+= delta*f(a+i*delta)
		count+=1
	end
	return sum, count
end

function newton_two(f, a, b, delta)
	N = (b-a)/(1.0*delta)
	h = delta/2.0
	sum = h*f(a)/3+h*f(b)/3
	count=2
	for n in 1:N-1
		sum+=(h/3)*(2*f(a+2*n*h)+4*f(a+(2*n-1)*h))
		count+=2
	end
	sum+= 4*h*f(a+(2*N-1)*h)/3
	count+=1
	return sum, count
end


function relative_error(numerical, true_sol)
	numerator = abs(numerical-true_sol)
	return numerator/(1.0*true_sol)
end

N = [10,100,1000,10000,100000,1000000,10000000]
A = 0.0
B = 10.0
delts = [(B-A)/(1.0*n) for n in N]

function test(x)
	return sin(x^2)/sqrt(x^2+1)
end


exact_sol = .47519858913634741151

F_rectangular(delta) = newton_zero(test, A, B, delta)
F_trapezoid(delta) = newton_one(test, A, B, delta)
F_simpsons(delta) = newton_two(test, A, B, delta)

rect_values_err = [relative_error(newton_zero(test,A,B,d)[1], exact_sol) for d in delts]
rect_func_calls = [newton_zero(test,A,B,delta)[2] for delta in delts]

trap_values_err = [relative_error(newton_one(test,A,B,delta)[1], exact_sol) for delta in delts]
trap_func_calls = [newton_one(test,A,B,delta)[2] for delta in delts]

simp_values_err = [relative_error(newton_two(test,A,B,d)[1],exact_sol) for d in delts]
simp_func_calls = [newton_two(test,A,B,delta)[2] for delta in delts]

rich_rect = [richardson_extrapolate(F_rectangular, d, 2, 1) for d in delts]
rich_rect_err = [relative_error(r[1], exact_sol) for r in rich_rect]
rich_rect_calls = [r[2] for r in rich_rect]

rich_trap = [richardson_extrapolate(F_trapezoid, d, 3, 2) for d in delts]
rich_trap_err = [relative_error(r[1], exact_sol) for r in rich_trap]
rich_trap_calls = [r[2] for r in rich_trap]

rich_simp = [richardson_extrapolate(F_simpsons, d, 2, 4) for d in delts]
rich_simp_err = [relative_error(r[1], exact_sol) for r in rich_simp]
rich_simp_calls = [r[2] for r in rich_simp]


using PyPlot
figure()
loglog(rect_func_calls,rect_values_err, label="Rectangular Rule")
loglog(rich_rect_calls,rich_rect_err,label="Richardson's Extrapolation Rectangular Rule")

legend()
xlabel("Number of Function Calls")
ylabel("Relative Error")
title("Problem 2 Rectangular Rule")
savefig("problem2_rect.png")

using PyPlot
figure()
loglog(trap_func_calls, trap_values_err,label="Trapezoid Rule")
loglog(rich_trap_calls,rich_trap_err,label="Richardson's Extrapolation Trapezoid Rule")
legend()
xlabel("Number of Function Calls")
ylabel("Relative Error")
title("Problem 2 Trapezoid Rule")
savefig("problem2_trap.png")

using PyPlot
figure()
loglog(simp_func_calls,simp_values_err, label="Simpson's Rule")
loglog(rich_simp_calls,rich_simp_err,label="Richardson's Extrapolation Simpson's Rule")
legend()
xlabel("Number of Function Calls")
ylabel("Relative Error")
title("Problem 2 Simpson's Rule")
savefig("problem2_simp.png")
println("done")