using LinearAlgebra

function inv_iter(A, b, mu; num_iter=10000, ender=sqrt(eps()))
    #make special A matrix
    special_A = inv(A-mu*Diagonal(ones(size(A)[1])))

    for i in 1:num_iter
        # calculate the matrix-by-vector product Ab
        b_next = special_A*b

        # calculate the norm
        b_next_norm = b_next/norm(b_next)

        # Termination criterion
        if norm(b - b_next_norm) < ender*norm(b_next_norm)
            break
        end 

        b = b_next_norm
    end
    return transpose(b)*A*b, b
end


function all_evecs(A)
    b0 = [1.0 2.0 3.0]
    mus = eigvals(A)-.000001*ones(size(A)[1])
    evecs = []
    evals = []

    for i in 1:size(A)[1]
        val, vec = inv_iter(A, transpose(b0), mus[i])
        push!(evecs, vec)
        push!(evals, val)
    end
    return evecs, evals
end



A = [2.0 -3.0 1.0;
    -3.0 1.0 4.0;
    2.0 4.0 -1.0]

evecs, evals = all_evecs(A)
println(evals)
println(evecs)

#Validate that they are eigenvectors
evec1=evecs[1]
evec2=evecs[2]
evec3=evecs[3]

println(evec1*evals[1])
println(A*evec1)

println(evec2*evals[2])
println(A*evec2)

println(evec3*evals[3])
println(A*evec3)


