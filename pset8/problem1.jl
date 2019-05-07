function DCT(f)
    N = length(f)
    fhat = zeros(ComplexF64, N)
    ghat = zeros(ComplexF64, 2*N)
    nus = zeros(ComplexF64,N)
    gnus = zeros(ComplexF64, 2*N)
    for nu=-N/2:N/2-1
        nus[Integer(N/2+nu+1)]=nu
        for n=0:N-1
            fhat[Integer(N/2+nu+1)] += f[n+1]*exp(-2im*pi*n*nu/N)/N #cos(2*pi*n*nu/N)/N
        end
    end
    for nu=-N:N-1
        gnus[Integer(N+nu+1)]=nu
        for n=0:N-1
            ghat[Integer(N+nu+1)] += (f[(n)+1]+f[(N-1-n)+1])*cos(pi*n*nu/N)/(2*N)
        end
    end    
    return fhat, ghat, nus, gnus
end


f = [0,1,2,3,4,5,6,7,8,9]
fh, gh, nus, gnus = DCT(f)
println(length(fh))
println(length(gh))
println(length(nus))

using PyPlot
figure()
plot(nus, fh, "o-")
title("Problem 1 - DFT for f")
xlabel("nu")
ylabel("hatf_nu")
savefig("problem1_f.png")
figure()
plot(gnus, gh, "o-")
title("Problem 1 - DFT for g")
xlabel("nu")
ylabel("hatg_nu")
savefig("problem1_g.png")


print("done")


