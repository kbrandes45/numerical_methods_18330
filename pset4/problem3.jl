using LinearAlgebra

function second_derivative(N, h, omega)
    # construct the matrix A implementing the fourth derivative

    A = diagm( 0 => (omega^2-2.0/(h^2.0))*ones(N),
               1 => (1/h^2.0)*ones(N-1),
              -1 => (1/h^2)*ones(N-1))
    return A
end



function forcing(t)
	return cos(t^2)
end

function solve(N, t_0, a, t_1, b, omega, forcing, makeA)
	h=(t_1-t_0)/(N+1)
	fdoubleprime = zeros(N)
	sparse_bc = zeros(N)
	sparse_bc[1]=a
	sparse_bc[N]=b
	time = zeros(N)
	for i in 1:N
		t = i*h
		fdoubleprime[i]=forcing(t)
		time[i] = t
	end
	q = fdoubleprime-(sparse_bc)/(h^2.0)
	A = makeA(N, h, omega)
	f = A \ q
	return f, time, fdoubleprime
end

N=1000
f, time, forcing_traj =  solve(N, 0, .3, 10, -2.9, 4.3, forcing, second_derivative)
println(f)


using PyPlot
figure()
plot(time, f , label="x(t) Trajectory")
plot(time, forcing_traj , label="Forcing Trajectory")

legend()
xlabel("Time t")
title("Problem 3 - Boundary Value Problem")
savefig("problem3_partb.png")
print("done")



function forcing2(t)
	return cos(4.3*t)
end

N=1000
f, time, forcing_traj =  solve(N, 0, .3, 10, -2.9, 4.3, forcing2, second_derivative)
println(f)


using PyPlot
figure()
plot(time, f , label="x(t) Trajectory")
plot(time, forcing_traj , label="Forcing Trajectory")

legend()
xlabel("Time t")
title("Problem 3 - Boundary Value Problem")
savefig("problem3_partc.png")
print("done")

