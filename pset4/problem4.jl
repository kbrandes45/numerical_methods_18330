using LinearAlgebra

function lx_matrix(N,h)
    A = diagm( 0 => (-2.0/h^2)*ones(N))
    C = diagm( 0 => (1.0/h^2)*ones(N))
    B = diagm(0=>zeros(N))
    final_matrix = zeros(N, N^2)
    for i in 1:N
        first_row = diagm(0=>zeros(N))
        for j in 1:N
            if (i==1)
                if (j==1)
                    first_row = A
                elseif (j==2)
                    first_row = [first_row C]
                else
                    first_row = [first_row B]
                end
            else
                if (j==1)
                    if (j==i-1)
                        first_row = C
                    elseif (j==i)
                        first_row = A
                    elseif (j==i+1)
                        first_row =  C
                    else
                        first_row = B
                    end
                else
                    if (j==i-1)
                        first_row = [first_row C]
                    elseif (j==i)
                        first_row = [first_row A]
                    elseif (j==i+1)
                        first_row = [first_row C]
                    else
                        first_row = [first_row B]
                    end
                end
            end
        end
        if (i==1)
            final_matrix = first_row
        else
            final_matrix = [final_matrix; first_row]
        end
    end
    return final_matrix
end



function ly_matrix(N, h)
    A = diagm( 0 => (-2.0/h^2)*ones(N),
               1 => (1/h^2)*ones(N-1),
              -1 => (1/h^2)*ones(N-1))
    B = diagm(0=>zeros(N))
    final_matrix = zeros(N, N^2)
    for j in 1:N
        first_row = zeros(N,N)
        for i in 1:N
            if (i==1)
                if (i==j)
                    first_row = A
                else
                    first_row = B
                end
            else
                if (i==j)
                    first_row = [first_row A]
                else
                    first_row = [first_row B]
                end
            end
        end
        if (j == 1)
            final_matrix = first_row
        else
            final_matrix = [final_matrix; first_row]
        end
    end
    return final_matrix
end


function solve(N, rho)
    h = 1/(N+1)
    rho_vals = zeros(N^2)
    for i in 1:N
        #x loop=x
        x_val = 0+i*h
        for j in 1:N
            #y loop = j
            y_val = 0+j*h
            index = N*(i-1)+j
            rho_vals[index] = rho(x_val, y_val)
        end
    end
    L = lx_matrix(N,h)+ly_matrix(N,h)
    phi = L\rho_vals
    return phi
end


function dx_matrix(N)
    A=diagm(0=>ones(N))
    B=diagm(0=>-1*ones(N))
    C=zeros(N,N)
    final_matrix = zeros(N, N^2)
    for i in 1:N
        first_row = diagm(0=>zeros(N))
        for j in 1:N
            if (j==1)
                if (j==i-1)
                    first_row = B
                elseif (j==i)
                    first_row = C
                elseif (j==i+1)
                    first_row =  A
                else
                    first_row = C
                end
            else
                if (j==i-1)
                    first_row = [first_row B]
                elseif (j==i)
                    first_row = [first_row C]
                elseif (j==i+1)
                    first_row = [first_row A]
                else
                    first_row = [first_row C]
                end
            end
            
        end
        if (i==1)
            final_matrix = first_row
        else
            final_matrix = [final_matrix; first_row]
        end
    end    
    return final_matrix
end

function dy_matrix(N)
    A=diagm(1=>ones(N-1),
            -1=>-1*ones(N-1))
    B=zeros(N,N)
    final_matrix = zeros(N, N^2)
    for j in 1:N
        first_row = zeros(N,N)
        for i in 1:N
            if (i==1)
                if (i==j)
                    first_row = A
                else
                    first_row = B
                end
            else
                if (i==j)
                    first_row = [first_row A]
                else
                    first_row = [first_row B]
                end
            end
        end
        if (j == 1)
            final_matrix = first_row
        else
            final_matrix = [final_matrix; first_row]
        end
    end
    return final_matrix
end



function electric(N, phi_vals)
    h=1/(N+1)
    D = (1/(2*h))*dx_matrix(N)+(1/(2*h))*dy_matrix(N)
    return -D*phi_vals
end

function rho(x,y)
    return exp(-50*(x-.5)^2-50*(y-.5)^2)
end

function rho2(x,y)
    return exp(-100*(x-.25)^2-100*(y-.25)^2)-exp(-100*(x-.75)^2-100*(y-.75)^2)
end

N=20
phi = solve(N, rho)
res = electric(N, phi)

phi2 = solve(N,rho2)
res2 = electric(N,phi2)

function convert(inp,N)
    final = reshape(inp, N, N)
    return final
end

phi_mat = convert(phi,N)
println(phi_mat)

using PyPlot
figure()
imshow(phi_mat)
colorbar()
xlabel("X")
ylabel("Y")
title("Electrostatic Potential for rho 1")
savefig("problem4_phi1.png")

rho_mat = convert(res,N)
println(rho_mat)

using PyPlot
figure()
imshow(rho_mat)
colorbar()
xlabel("X")
ylabel("Y")
title("Electric Field for rho 1")
savefig("problem4_rho1.png")


phi_mat = convert(phi2,N)
println(phi_mat)

using PyPlot
figure()
imshow(phi_mat)
colorbar()
xlabel("X")
ylabel("Y")
title("Electrostatic Potential for rho 1")
savefig("problem4_phi2.png")



rho_mat = convert(res2,N)
println(rho_mat)

using PyPlot
figure()
imshow(rho_mat)
colorbar()
xlabel("X")
ylabel("Y")
title("Electric Field for rho 1")
savefig("problem4_rho2.png")

println("done")


