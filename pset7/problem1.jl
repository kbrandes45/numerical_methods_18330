using LinearAlgebra

function integrate_mc(f, a, b, N)
    D = length(a)
    println(D)
    S = 0.0
    
    # volume of the hypercube V = (b_1 - a_1)×(b_2 - a_2)× ... × (b_D - a_D)
    V = prod(b .- a)

    for i=1:N
        # sample a point uniformly from [0, 1]×...×[0,1]
        x = rand(D)

        # scale to [a[1], b[1]] × ... × [a[D], b[D]]
        x = a .+ (b .- a).*x
        
        S += f(x)
    end
    
    return S.* (V/N)
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
samples=10000

res = integrate_mc(f,a,b,samples)
print(res)


samples = [10, 100, 1000, 10000, 100000, 1000000, 10000000]
res = [integrate_mc(f,a,b,s) for s in samples]

function num_digits_correct(x)
    println(x)
    println(pi/4.0)
    err = abs(x-pi/4.0)
    val = log10(err)
    println(val)
    digits = round(-1.0*val)
    println(digits)
    return digits
end
correct_digits = [num_digits_correct(x) for x in res]


using PyPlot
figure()
semilogx(samples,correct_digits)

xlabel("Number of Samples")
ylabel("Number of Correct Digits")
title("Problem 1 Estimate Pi")
savefig("problem1.png")
print("done")