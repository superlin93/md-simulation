function cart_pos = frac_to_cart(frac_pos, lattice)
% Calculate cartesian coordinates from the fractional position
    cart_pos = zeros(3,1);
    for i = 1:3
        cart_pos(i) = lattice(1,i)*frac_pos(1) + lattice(2,i)*frac_pos(2) + lattice(3,i)*frac_pos(3); 
    end
end