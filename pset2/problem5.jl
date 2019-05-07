function planets(planetx, planety, x1, y1, x2, y2)
	dist1=((planetx-x1)^2+(planety-y1)^2)^(.5)
	dist2=((planetx-x2)^2+(planety-y2)^2)^(.5)
	term1x = -1.0*(planetx-x1)/(dist1^3)+-1.0*(planetx-x2)/(dist2^3)
	term1y = -1.0*(planety-y1)/(dist1^3)+-1.0*(planety-y2)/(dist2^3)
	return [term1x,term1y]
end


function leapfrog(starting_x1, starting_x1_prime, starting_x2, starting_x2_prime, starting_x3, starting_x3_prime, h, N, f)

	planet_1_x_list = [starting_x1]
	planet_1_y_list = [starting_x1_prime]
	planet_2_x_list = [starting_x2]
	planet_2_y_list = [starting_x2_prime]
	planet_3_x_list = [starting_x3]
	planet_3_y_list = [starting_x3_prime]
	for i in 1:N
		f_x_list = f(planet_1_x_list[i][1], planet_1_x_list[i][2], planet_2_x_list[i][1], planet_2_x_list[i][2], planet_3_x_list[i][1], planet_3_x_list[i][2])
		v_half_p1x = planet_1_y_list[i][1]+h/2.0*f_x_list[1]
		v_half_p1y = planet_1_y_list[i][2]+h/2.0*f_x_list[2]

		x_next_p1x = planet_1_x_list[i][1]+h*v_half_p1x
		x_next_p1y = planet_1_x_list[i][2]+h*v_half_p1y

		f_x_list2 = f(planet_2_x_list[i][1], planet_2_x_list[i][2], planet_1_x_list[i][1], planet_1_x_list[i][2], planet_3_x_list[i][1], planet_3_x_list[i][2])
		v_half_p2x = planet_2_y_list[i][1]+h/2.0*f_x_list2[1]
		v_half_p2y = planet_2_y_list[i][2]+h/2.0*f_x_list2[2]

		x_next_p2x = planet_2_x_list[i][1]+h*v_half_p2x
		x_next_p2y = planet_2_x_list[i][2]+h*v_half_p2y		

		f_x_list3 = f(planet_3_x_list[i][1], planet_3_x_list[i][2], planet_2_x_list[i][1], planet_2_x_list[i][2], planet_1_x_list[i][1], planet_1_x_list[i][2])
		v_half_p3x = planet_3_y_list[i][1]+h/2.0*f_x_list3[1]
		v_half_p3y = planet_3_y_list[i][2]+h/2.0*f_x_list3[2]

		x_next_p3x = planet_3_x_list[i][1]+h*v_half_p3x
		x_next_p3y = planet_3_x_list[i][2]+h*v_half_p3y

		f_xnew_list = f(x_next_p1x, x_next_p1y, x_next_p2x, x_next_p2y, x_next_p3x, x_next_p3y)
		y_next_p1x = v_half_p1x+h/2.0*f_xnew_list[1]
		y_next_p1y = v_half_p1y+h/2.0*f_xnew_list[2]

		f_xnew_list2 = f(x_next_p2x, x_next_p2y, x_next_p1x, x_next_p1y, x_next_p3x, x_next_p3y)
		y_next_p2x = v_half_p2x+h/2.0*f_xnew_list2[1]
		y_next_p2y = v_half_p2y+h/2.0*f_xnew_list2[2]

		f_xnew_list3 = f(x_next_p3x, x_next_p3y, x_next_p2x, x_next_p2y, x_next_p1x, x_next_p1y)
		y_next_p3x = v_half_p3x+h/2.0*f_xnew_list3[1]
		y_next_p3y = v_half_p3y+h/2.0*f_xnew_list3[2]

		push!(planet_1_x_list, [x_next_p1x, x_next_p1y])
		push!(planet_1_y_list, [y_next_p1x, y_next_p1y])

		push!(planet_2_x_list, [x_next_p2x, x_next_p2y])
		push!(planet_2_y_list, [y_next_p2x, y_next_p2y])

		push!(planet_3_x_list, [x_next_p3x, x_next_p3y])
		push!(planet_3_y_list, [y_next_p3x, y_next_p3y])

	end
	return planet_1_x_list, planet_2_x_list, planet_3_x_list
end

p1, p2, p3 = leapfrog([-.7,.35], [.99,.078],[1.1,-.07],[.1,.47],[-.4,-.3],[-1.1,-.53],.1,20,planets)
#p1, p2, p3 = leapfrog([-.5,.35], [.99,.078],[1.35,-.07],[.1,.47],[-.5,-.35],[-1.1,-.53],.1,20,planets)
#p1, p2, p3 = leapfrog([-.7,.35], [.99,.108],[1.1,-.07],[.125,.4],[-.4,-.3],[-1.4,-.61],.1,20,planets)
planet_1_x = [coord[1] for coord in p1]
planet_1_y = [coord[2] for coord in p1]

planet_2_x = [coord[1] for coord in p2]
planet_2_y = [coord[2] for coord in p2]

planet_3_x = [coord[1] for coord in p3]
planet_3_y = [coord[2] for coord in p3]


using PyPlot
figure()
plot(planet_1_x, planet_1_y, "o-", label=L"Planet 1 Motion")
plot(planet_2_x, planet_2_y,"o-",label=L"Planet 2 Motion")
plot(planet_3_x,planet_3_y,"o-",label=L"Planet 3 Motion")
xlabel("X Values")
ylabel("Y Values")
legend()
title("Planetary Motion")
savefig("problem5planets4.png")
print("planets done")


