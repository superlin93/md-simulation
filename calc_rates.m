function [jump_names, jumps_total, rates, e_act] = calc_rates(sites, sim_data)
% Calculate the jump rates between sites (and the standard deviation) and activation energies
    nr_stable = size(sites.succes,1);
    nr_diff_names = size(sites.names,1);
% In how many parts the simulation is separated for calculating the standard deviation:
    nr_parts = sites.nr_parts;

%% Find the stable names, transition states start with 'Transition'
    nr_stable_names = 0;
    for i = 1:nr_diff_names
       test = strsplit(sites.names{i}, ':');
       if ~strcmp('Transition', test{1})    
            nr_stable_names = nr_stable_names + 1;            
       end
    end 
    nr_combi = nr_stable_names^2;   
    
%% The possible name-combinations
    % allow for A-B and B-A to check if they are equal 
    counter = 0;
    jump_names = cell(nr_combi,1); 
    for site1 = 1:nr_stable_names %loop over stable names, are at the start of grouped.names
        for site2 = 1:nr_stable_names  
            counter = counter + 1;
            jump_names{counter} = [sites.names{site1}, '_to_', sites.names{site2}];
        end
    end  
    
   % Find the succesfull jumps and combine them by name
    jumps = zeros(nr_combi, nr_parts);
    for site1 = 1:nr_stable %loop over succes
        for site2 = 1:nr_stable
            for part = 1:nr_parts
                if sites.succes(site1, site2, part) > 0
                % Get the name of this combination    
                    name = [sites.site_names{site1}, '_to_', sites.site_names{site2}];
                % Find to what this corresponds to in jumps()
                    name_nr = 0;
                    for i = 1:nr_combi
                        if strcmp(jump_names{i}, name)
                            name_nr = i;
                        end
                    end
                % Add it to the correct name:
                    if name_nr ~= 0 % QUICK FIX, should not be necessary!!
                        jumps(name_nr, part) = jumps(name_nr, part) + sites.succes(site1, site2, part);
                    else
                        disp('ERROR! No name found for a jump! (in calc_rates.m)')
                    end
                end
            end
        end
    end
    
    %% Calculate activation energies, jumps rates per atom per second, and the standard deviation:
    rates = zeros(nr_combi,2);
    e_act = zeros(nr_combi,1);
    nr_diff_atoms = size(sites.atoms,2); %nr of diffusing atoms
    jump_freqs = jumps./(nr_diff_atoms*(sim_data.total_time/nr_parts)); %obtain the jump frequency per diffusing atom
    jumps_total = sum(jumps, 2);  
    for i = 1:nr_combi
        % The average
        rates(i,1) = mean(jump_freqs(i,:));
        % The standard deviation
        rates(i,2) = std(jump_freqs(i,:));
    end
    
    % The activation energy in eV (ASSUMING IONS WITH A +1 CHARGE NOW!)
    for i = 1:nr_combi
        if rates(i,1) > 0
            % Split jump_names        
            names = strsplit(jump_names{i}, '_to_');
            for j = 1:nr_stable_names % Find the occupancy of the starting site:
                if strcmp(sites.stable_names{j}, names{1})
                    atom_percentage = sites.atom_locations(j);
                end
            end
     % The effective rate, i.e. how long an atom on average stays at this site:
            eff_rate = jumps_total(i)/(atom_percentage*sim_data.total_time*sim_data.nr_diffusing);
            if strcmp(names{1}, names{2})
            % For A-A jumps divide by two for a fair comparison of A-A jumps vs. A-B and B-A)
                eff_rate = eff_rate/2; 
            end
            % The average activation energy:
            e_act(i,1) = -log(eff_rate/sim_data.attempt_freq) * ...
                (sim_data.k_boltzmann * sim_data.temperature)/sim_data.e_charge; 
        else %When no jumps occur the activation energy is unknown:
            e_act(i,:) = NaN;
        end
    end
    
end