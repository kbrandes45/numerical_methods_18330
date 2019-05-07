println("start")
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

function precent_different(val1, val2)
	num = abs(val1-val2)
	denom = (val1+val2)/2.0
	return (num/denom)*100
end

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
		diff = precent_different(ie_step, rk4_step)
		print("diff ")
		println(diff)

		if (diff <= 1.0)
			push!(h_vals, h_i*1.2)
			push!(u_vals, rk4_step)
			push!(time, t_i+h_i)
		elseif (diff >= 15.0)
			h_vals[length(h_vals)]= h_i*.8
		else
			push!(h_vals, h_i)
			push!(u_vals, rk4_step)
			push!(time, t_i+h_i)
		end
	end
	
	return time, u_vals, h_vals
end	


function test(t,u)
	return cos(t*u^2)
end

function test2(t,u)
	return exp(t)
end

time, u_vals, h_vals = adaptive_ode_solver(test, 0.0, 50, 3.0, .5)
using PyPlot
figure()
semilogy(time, u_vals , label="Solution u(t) curve")
semilogy(time, h_vals,  label="H step sizes")
legend()
xlabel("Time")
title("Problem 1 - Adaptive ODE Solver")
savefig("problem1_ode.png")
print("done")


