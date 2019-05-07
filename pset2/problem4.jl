function x(last_x, h, last_x_prime)
	return last_x+h*last_x_prime
end

function xprime(last_x, h, last_x_prime)
	return last_x_prime-h*last_x
end

function euler(starting_x, starting_x_prime, h, N)
	x_list = [starting_x*1.0]
	y_list = [starting_x_prime*1.0]
	for i in 1:N
		x_next = x(x_list[i], h, y_list[i])
		x_prime_next = xprime(x_list[i], h, y_list[i])
		push!(x_list, x_next)
		push!(y_list, x_prime_next)
	end
	return x_list, y_list
end

function exact(starting_x, starting_x_prime, h, N)
	#exact soln x(t)=sin(x)
	x_list = [starting_x*1.0]
	y_list = [starting_x_prime*1.0]
	t=0
	for i in 1:N
		x_next = sin(t)
		y_next = cos(t)
		t+=h
		push!(x_list, x_next)
		push!(y_list,y_next)
	end
	return x_list, y_list
end

function energy(x_list, y_list, h, N)
	time = []
	energy = []
	t = 0.0
	for i in 1:N
		push!(time, t)
		e = .5*x_list[i]^2+.5*y_list[i]^2
		push!(energy, e)
		t += h
	end
	return time, energy

end

#plot x vs x_prime
x_list, y_list = euler(0, 1, .05, 20)
x_list_exact, y_list_exact = exact(0,1,.05,20)
using PyPlot
figure()
plot(x_list,y_list,"o-",label=L"Euler Method")
plot(x_list_exact, y_list_exact,"o-",label=L"Exact Solution")
xlabel("X values")
ylabel("X' values")
legend()
title("Problem 2 - Second Order ODE")
savefig("problem4.png")
print("half done")

x_list, y_list = euler(0, 1, .2, 200)
time, en = energy(x_list, y_list, .2, 200)

using PyPlot
figure()
plot(time,en)
xlabel("Time")
ylabel("Energy")
title("Euler Method for Conservation of Energy")
savefig("problem4energy.png")
print("done")

function leapfrog(starting_x, starting_x_prime, h, N)
	x_list = [1.0*starting_x]
	y_list = [1.0*starting_x_prime]
	for i in 1:N
		v_half = y_list[i]+h/2.0*-(x_list[i])
		x_next = x_list[i]+h*v_half
		y_next = v_half+h/2.0*(-x_next)
		push!(x_list, x_next)
		push!(y_list, y_next)
	end
	return x_list, y_list
end

x_list, y_list = leapfrog(0,1.0, .05, 20)
x_list1, y_list1 = leapfrog(0,1.0, .2, 200)

time, en = energy(x_list1, y_list1, .2, 200)
x_list_exact, y_list_exact = exact(0,1,.05,20)
using PyPlot
figure()
plot(x_list, y_list, "o-", label=L"Leapfrog Method")
plot(x_list_exact, y_list_exact,"o-",label=L"Exact Solution")
xlabel("X Values")
ylabel("X' Values")
legend()
title("Leapfrog Solution to Second Order ODE")
savefig("problem4leapfrog.png")
print("leapfrog done")

using PyPlot
figure()
plot(time, en)
xlabel("Time")
ylabel("Energy")
title("Leapfrog Solution Energy")
savefig("problem4leapfrogenergy.png")
print("really done")


#Forward vs reverse integration
x_list, y_list = euler(0, 1, .05, 20)
print(x_list[20],y_list[20])
x_rev, y_rev = euler(0.8417903781742831, 0.5402146250460914,-.05,20)
using PyPlot
figure()
plot(x_list,y_list,"o-",label=L"Forward Euler")
plot(x_rev, y_rev,"o-",label=L"Reversed Euler")
plot(x_list_exact, y_list_exact,"o-",label=L"Exact Solution")
xlabel("X values")
ylabel("X' values")
legend()
title("Problem 4 - Forward vs Reverse Euler")
savefig("problem4forwardreverse.png")




#Forward vs reverse integration
x_list, y_list = leapfrog(0, 1, .05, 20)
print(x_list[21],y_list[21])
x_rev, y_rev = leapfrog(0.8417903781742831, 0.5402146250460914,-.05,20)
using PyPlot
figure()
plot(x_list,y_list,"o-",label=L"Forward Leapfrog")
plot(x_rev, y_rev,"o-",label=L"Reversed Leapfrog")
plot(x_list_exact, y_list_exact,"o-",label=L"Exact Solution")
xlabel("X values")
ylabel("X' values")
legend()
title("Problem 4 - Forward vs Reverse Leapfrog")
savefig("problem4leapfrog_frev.png")
