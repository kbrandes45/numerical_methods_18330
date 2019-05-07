
using LinearAlgebra
function newton_nd(f, J, x0; N=10)
    x = x0
    for i=1:N
        x = x .- (J(x) \ f(x))
    end 
    return x
end

function get_tolerance(newtonND, f, J, x0)
    n=1
    x_N=newtonND(f,J,x0; N=1)
    while norm(f(x_N))>10^(-8.0)
    	n+=1
    	x_N=newtonND(f,J,x0;N=n)
    end
    return x_N, n, norm(f(x_N))
end




#Part A
x0 = [1.0; 1.0] #x1,x2
f(x) = [3*x[1]^2 - x[2]^2; 3*x[1]*x[2]^2-x[1]^3-1]
println("part a")
J(x) = [
    (6*x[1])   		(-2*x[2])
    (3*x[2]^2-3*x[1]^2) (6*x[1]*x[2])
]


#x_N = newton_nd(f, J, x0; N=10)
x_N, iters, tol = get_tolerance(newton_nd, f, J, x0)
println(x_N)
println(iters)
println(tol)



#Part B
x0 = [1.0; 1.0; 1.0] #x1,x2
f(x) = [6*x[1] - 2*cos(x[2]*x[3])-1; 9*x[2]+(x[1]^2+sin(x[3])+1.06)^(.5)+.9; 60*x[3]+3*exp(-x[1]*x[2])+10*pi-3]
println("part b")
J(x) = [
    (6) (2*x[3]*sin(x[2]*x[3])) (2*x[2]*sin(x[2]*x[3]))
    ((x[1])/((x[1]^2+sin(x[3])+1.06)^(.5))) (9) (cos(x[3])/(2*(x[1]^2+sin(x[3])+1.06)^(.5)))
    (-3*x[2]*exp(-x[1]*x[2])) (-3*x[1]*exp(-x[1]*x[2])) (60)
]


x_N, iters, tol = get_tolerance(newton_nd, f, J, x0)
println(x_N)
println(iters)
println(tol)

