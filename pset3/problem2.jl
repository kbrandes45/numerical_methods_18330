print("starting")

function xfunc(t)
	return sin(t)
end

function xprime(t)
	return cos(t)
end

x = [xfunc(t) for t in 1:10]
y = [xprime(t) for t in 1:10]

using PyPlot
figure()
plot(x,y)
xlabel("X")
ylabel("X'")
title("Problem 2 - Analytical Solution")
savefig("problem_1_analytical.png")

function improved_euler_onestep(f, t_n, u_n, h)
    u_next_euler = u_n + h*f(t_n+h, u_n)
    u_next = u_n + h/2.0*(f(t_n+h, u_n) + f(t_n+h, u_next_euler))
    return u_next
end

function rk4_onestep(f, t_n, u_n, h)
	k1 = h*f(t_n, u_n)
	k2 = h*f(t_n+h/2.0, u_n+k1/2.0)
	k3 = h*f(t_n+h/2.0, u_n+k2/2.0)
	k4 = h*f(t_n+h, u_n+k3)
	return u_n+1/6.0*(k1+2*k2+2*k3+k4)
end

function precent_different_vec(val1, val2)
	num = (abs(val1[1]-val2[1])+abs(val1[2]-val2[2]))/2.0
	denom = ((val1[1]+val2[1])/2.0+(val1[2]+val2[2])/2.0)/2.0
	return (num/denom)*100
end

function precent_diff_vec(val1, val2)
	vec_diff1 = val1[1]-val2[1]
	vec_diff2 = val1[2]-val2[2]
	num = (vec_diff1^2+vec_diff2^2)^(1/2.0)
	avg_vec_diff1 = (val1[1]+val2[1])/2.0
	avg_vec_diff2 = (val1[2]+val2[2])/2.0
	denom = (avg_vec_diff1^2+avg_vec_diff2^2)^(1/2.0)
	return (num/denom)*100
end



println(precent_diff_vec([1.0,1.0],[1.0,1.0]))
println(precent_diff_vec([1.0,2.0],[1.0,1.0]))

println(precent_different_vec([1.0,2.0],[1.0,1.0]))






function adaptive_ode_solver(f, t_0, t_f, u_0, h)
	time = [t_0]
	u_vals = [u_0]
	h_vals = [h]

	while time[length(time)]<t_f
		t_i = time[length(time)]
		u_i = u_vals[length(u_vals)]
		h_i = h_vals[length(h_vals)]
		println(h_i)
		print("time ")
		println(t_i)

		rk4_step = rk4_onestep(f, t_i, u_i, h_i)
		print("rk4_step ")
		println(rk4_step)
		ie_step = improved_euler_onestep(f, t_i, u_i, h_i)
		print("ie_step ")
		println(ie_step)
		diff = precent_diff_vec(ie_step, rk4_step)
		println(diff)

		if (diff <= 1.0)
			push!(h_vals, h_i*1.2)
			push!(u_vals, rk4_step)
			push!(time, t_i+h_i)
		elseif (diff >= 25.0)
			h_vals[length(h_vals)]= h_i*.8
		else
			push!(h_vals, h_i)
			push!(u_vals, rk4_step)
			push!(time, t_i+h_i)
		end
	end
	
	return time, u_vals, h_vals
end	

#u is a vector, t is a scalar
function simplified_func(t, u)
	x_1 = u[2]
	x_2 = 5*(1-u[1]^2)*u[2]-u[1]
	return [x_1, x_2]
end

function parta_func(t, u)
	x_1 = u[2]
	x_2 = 0*(1-u[1]^2)*u[2]-u[1]
	return [x_1, x_2]
end

function nonzero_f(t)
	return pi/2.0+atan(1000*sin(t))
end

function partc_func(t,u)
	x_1 = u[2]
	x_2 = 5*(1-u[1]^2)*u[2]-u[1]+nonzero_f(t)
	return [x_1, x_2]
end

#println(improved_euler_onestep(simplified_func,0.0,[0.0,1.0],.05))

#println(rk4_onestep(simplified_func,0.0,[0.0,1.0],.05))

time, u_vals, h_vals = adaptive_ode_solver(parta_func, 0.0, 40.0, [0.0,1.0], .5)
using PyPlot
figure()

x_vals = [u[1] for u in u_vals]
x_prime_vals = [u[2] for u in u_vals]

plot(time,x_vals , label="Solution x(t) curve")
plot(time, x_prime_vals , label="Solution x'(t) curve")
semilogy(time, h_vals,  label="H step sizes")
legend()
xlabel("Time")
title("Problem 2 - Adaptive ODE Solver with mu=0 and f(t)=0")
savefig("problem2_ode_a.png")


time, u_vals, h_vals = adaptive_ode_solver(simplified_func, 0.0, 40.0, [0.0,1.0], .5)
using PyPlot
figure()

x_vals = [u[1] for u in u_vals]
x_prime_vals = [u[2] for u in u_vals]

plot(time,x_vals , label="Solution x(t) curve")
plot(time, x_prime_vals , label="Solution x'(t) curve")
semilogy(time, h_vals,  label="H step sizes")
legend()
xlabel("Time")
title("Problem 2 - Adaptive ODE Solver with mu=5 and f(t)=0")
savefig("problem2_ode_b.png")

using PyPlot
figure()
plot(x_vals, x_prime_vals)
legend()
xlabel("X")
ylabel("X'")
title("Problem 2 - Phase Space with mu=5 and f(t)=0")
savefig("problem2_b_phase_space.png")


time, u_vals, h_vals = adaptive_ode_solver(partc_func, 0.0, 40.0, [0.0,1.0], 1.0)
using PyPlot
figure()

x_vals = [u[1] for u in u_vals]
x_prime_vals = [u[2] for u in u_vals]

plot(time,x_vals , label="Solution x(t) curve")
plot(time, x_prime_vals , label="Solution x'(t) curve")
semilogy(time, h_vals,  label="H step sizes")
legend()
xlabel("Time")
title("Problem 2 - Adaptive ODE Solver with mu=5 and f(t)=pi/2+arctan(1000*sin(t))")
savefig("problem2_ode_c.png")
print("done")



