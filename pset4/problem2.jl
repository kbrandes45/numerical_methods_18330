function f_func(x)
	return exp(x)
end

function g_func(x)
	return x^2
end

true_f = exp(1.0)
true_g = 2000.0


function relative_error(numerical, true_sol)
	numerator = abs(numerical-true_sol)
	return numerator/(1.0*true_sol)
end


function forward_diff_stencil(func, x, h)
	numerator = func(x+h)-func(x)
	return numerator/(1.0*h)
end

h_list = [1, .1, .01, .001, .0001, .00001, .000001, .0000001, .00000001, .000000001, .0000000001]

f_rel_err = [relative_error(forward_diff_stencil(f_func, 1.0, h), true_f) for h in h_list]
g_rel_err = [relative_error(forward_diff_stencil(g_func, 1000.0, h), true_g) for h in h_list]

using PyPlot
figure()
loglog(h_list, f_rel_err , label="Function f(x)=exp(x) at x=1.0")
loglog(h_list, g_rel_err,  label="Function g(x)=x^2 at x=1000.0")
legend()
xlabel("h")
ylabel("Relative Error")
title("Problem 3 - Relative Error for Forward Difference Stencil")
savefig("problem2_relerr.png")
print("done")



