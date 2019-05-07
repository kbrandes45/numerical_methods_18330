using LinearAlgebra
function eig_qr(A; N=1000, ϵ=√eps())
    """ Compute eigenvalues of A using the QR algorithm
    """
    F = qr(A)
    Q, R = F.Q, F.R
    eigenvec_mat = F.Q
    for i=1:N
        # multiply in "wrong" order
        A_new = R*Q

        # termination criterion
        if norm(A_new - A) < ϵ*norm(A_new)
            break
        end
        
        A = A_new
        
        # QR decomposition
        F = qr(A)
        Q, R = F.Q, F.R

        #increment eigenvector matrix
        eigenvec_mat = eigenvec_mat*Q

    end
    
    #get the eigenvalues from diag of A
    evals = zeros(size(A)[1])
    for i in 1:size(A)[1]
        evals[i] = A[i+(i-1)*size(A)[1]]
    end
    println("A matrix")
    println(A)
    return evals, eigenvec_mat
end

A = [2.0 -3.0 2.0;
    -3.0 1.0 4.0;
    2.0 4.0 -1.0]

evals, evecs = eig_qr(A)
println(evals)
println(evecs)

#Validate that they are eigenvectors
evec1=[evecs[1] evecs[2] evecs[3]]
evec2=[evecs[4] evecs[5] evecs[6]]
evec3=[evecs[7] evecs[8] evecs[9]]

println(evec1*evals[1])
println(A*transpose(evec1))

println(evec2*evals[2])
println(A*transpose(evec2))

println(evec3*evals[3])
println(A*transpose(evec3))

#Validate they are orthonormal
println(evec1*transpose(evec1)) #1
println(evec1*transpose(evec2)) #0
println(evec1*transpose(evec3)) #0

println(evec2*transpose(evec1)) #0
println(evec2*transpose(evec2)) #1
println(evec2*transpose(evec3)) #0

println(evec3*transpose(evec1)) #0
println(evec3*transpose(evec2)) #0
println(evec3*transpose(evec3)) #1



A = [2.0 -3.0 1.0;
    -3.0 1.0 4.0;
    2.0 4.0 -1.0]

evals, evecs = eig_qr(A)
println(evals)
println(evecs)

#Validate that they are eigenvectors
evec1=[evecs[1] evecs[2] evecs[3]]
evec2=[evecs[4] evecs[5] evecs[6]]
evec3=[evecs[7] evecs[8] evecs[9]]

println(evec1*evals[1])
println(A*transpose(evec1))

println(evec2*evals[2])
println(A*transpose(evec2))

println(evec3*evals[3])
println(A*transpose(evec3))

#Validate they are orthonormal
println(evec1*transpose(evec1)) #1
println(evec1*transpose(evec2)) #0
println(evec1*transpose(evec3)) #0

println(evec2*transpose(evec1)) #0
println(evec2*transpose(evec2)) #1
println(evec2*transpose(evec3)) #0

println(evec3*transpose(evec1)) #0
println(evec3*transpose(evec2)) #0
println(evec3*transpose(evec3)) #1