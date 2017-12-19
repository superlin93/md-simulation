function dist_sqrt = calc_dist_sqrd_frac(frac1, frac2, lattice)
% Given two fractional coordinates and the lattice calculate dist^2
% (^2 since taking the root takes a long time and is not always necessary)
    frac = frac1 - frac2;
% Perodic Boundary Conditions: 
    for i = 1:3
        if frac(i) > 0.5
            frac(i) = frac(i) - 1.0;
        elseif frac(i) < -0.5
            frac(i) = frac(i) + 1.0;
        end
    end
%  *lattice to get cartesian:
    cart = frac'*lattice;
% calculate distance^2
    dist_sqrt = cart(1)^2 + cart(2)^2 + cart(3)^2; 

end