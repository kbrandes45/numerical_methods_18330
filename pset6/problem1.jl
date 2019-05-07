function get_index(i,j,s)
	return i+s*(j-1)
end

function get_u(A,U,L,i,j,s)
	second_term = 0
	for k in 1:i-1
		second_term+=L[get_index(i,k,s)]*U[get_index(k,j,s)]
	end	

	return A[get_index(i,j,s)]-second_term
end

function get_l(A,U,L,i,j,s)
	second_term = 0
	for k in 1:j-1
		second_term+=L[get_index(i,k,s)]*U[get_index(k,j,s)]
	end	
	x= (A[get_index(i,j,s)]-second_term)/U[get_index(j,j,s)]
	return x
end


function lu_decomp(A)
	s = size(A)[1] 
	L = zeros(Float64, s, s)
	U = zeros(Float64, s, s)

	for i in 1:s
		#define u
		for j in 1:s
			if (i<=j)
				U[get_index(i,j,s)]=get_u(A,U,L,i,j,s)
			end
		end
		#define l
		for subi in 1:s
			if (i==subi)
				L[get_index(subi,i,s)]=1
			elseif (subi>i)
				L[get_index(subi,i,s)]=get_l(A,U,L,subi,i,s)
			end
		end
	end
	return L,U
end



function cholesky_decomp(A)
	s = size(A)[1] 
	L = zeros(Float64, s, s)

	#diagonals
	for i in 1:s
		for j in 1:s
			if (j==i)
				sums = 0
				for k in 1:i-1
					sums+=L[get_index(i,k,s)]^2
				end
				L[get_index(i,i,s)]=sqrt(A[get_index(i,i,s)]-sums)
			elseif (i<j)
				second_term = 0
				for k in 1:i-1
					second_term+=L[get_index(i,k,s)]*L[get_index(k,j,s)]
				end	
				L[get_index(j,i,s)]=(A[get_index(i,j,s)]-second_term)/L[get_index(i,i,s)]
			end
		end
	end
	return L, transpose(L)
end


function forward_sub(L,U,b)
	s = size(L)[1] 
	y = zeros(Float64,s)
	#solve l first 
	for i in 1:s
		sums = 0
		for k in 1:i-1
			sums += L[get_index(i,k,s)]*y[k]
		end
		y[i] = b[i]-sums
	end

	x = zeros(Float64,s)
	#solve u second for x vector via back sub
	for i in 1:s
		sums = 0
		j = s-i+1
		for k in j+1:s
			sums+=U[get_index(j,k,s)]*x[k]
		end
		x[j] = (y[j]-sums)/U[get_index(j,j,s)]
	end
	return x
	
end

using LinearAlgebra
println("Part A")
a_mat = [1.0 3.0 2.0;
		 2.4 -3.3 1.1;
		 -1.0 0.0 2.0]
out = lu_decomp(a_mat)
print(out[1])
print(out[2])
a_vec = [1.0 2.0 3.0]
println(forward_sub(out[1],out[2], a_vec))
#println(a_mat\transpose(a_vec))

println("Part B")
b_mat = [5.0 -3.0 2.0;
		-3.0 6.0 -1.0;
		2.0 -1.0 5.0]
b_vec = [1.0 2.0 3.0]
out = cholesky_decomp(b_mat)
println(out[1])
println(forward_sub(out[1],out[2], b_vec))
#println(cholesky(b_mat))

println("Part C")
c_mat = [2.0 -3.0 2.0;
		-3.0 3.0 -1.0;
		2.0 -1.0 2.0]
c_vec = [1.0 2.0 3.0]
#out = cholesky_decomp(c_mat)
out = lu_decomp(c_mat)
print(out[1])
print(out[2])
println(forward_sub(out[1],out[2], c_vec))


