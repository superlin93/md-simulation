%%%%%%%%%%%%%%%%%%%%%%
function sim_data = frac_and_disp(sim_data)
% Calculate cartesian coordinates and the displacement of the atoms...
    disp('Transforming cartesian to fractional coordinates and calculating displacement')
    sim_data.frac_pos = zeros(3, sim_data.nr_atoms, sim_data.nr_steps);
    sim_data.displacement = zeros(sim_data.nr_atoms, sim_data.nr_steps);

    d = zeros(3,1);
    for atom = 1:sim_data.nr_atoms
        uc = [0 0 0]; % Keep track of the amount of unit cells the atom has moved
        %For the first time step:
        sim_data.frac_pos(:,atom,1) = sim_data.cart_pos(:,atom,1)'/sim_data.lattice;
        start = sim_data.cart_pos(:,atom,1); % Starting position
        for time = 2:sim_data.nr_steps
            sim_data.frac_pos(:,atom,time) = sim_data.cart_pos(:,atom,time)'/sim_data.lattice;  
            for i = 1:3
                frac_diff = sim_data.frac_pos(i,atom,time)-sim_data.frac_pos(i,atom,time-1);
              % If frac_pos differs by more than 0.5 from the previous time step the atom changed unit cell
                if frac_diff > 0.5
                    uc(i) = uc(i)-1;
                elseif frac_diff < -0.5
                    uc(i) = uc(i)+1; 
                end
            end
            % Calculate displacement: 
            for i = 1:3
                d(i) = ( sim_data.cart_pos(i,atom,time) - start(i) + uc*sim_data.lattice(:,i) )^2;
            end
            sim_data.displacement(atom,time) = sqrt(sum(d));
        end
    end    

end