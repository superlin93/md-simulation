function [stable_names, sites_occup, atom_locations, occup_parts, atom_loc_parts] = calc_site_occups(sites)
% Determines what fraction of the time each type of site is occupied and the fraction of time atoms spent at a type of site
    nr_stable = size(sites.succes,1);
    nr_diff_names = size(sites.names,1);
    nr_steps = size(sites.atoms,1);
    nr_atoms = size(sites.atoms,2);
    
%% Find the stable names, transition states start with 'Transition'
    nr_stable_names = 0;
    for i = 1:nr_diff_names
       test = strsplit(sites.names{i}, ':');
       if ~strcmp('Transition', test{1})    
            nr_stable_names = nr_stable_names + 1; 
            stable_names(nr_stable_names,1) = test;
       end
    end 

    site_count = zeros(nr_stable_names,1); %How much of each site there are
    sites_occup = zeros(nr_stable_names,1); %How much all the sites are occupied (%)
    occup_parts = zeros(nr_stable, sites.nr_parts); %How much all the sites are occupied (%), at the different parts of the simulation.
    atom_locations = zeros(nr_stable_names,1); %What percentage of the time atoms are at a site
    for i = 1:nr_stable
        for j = 1:nr_stable_names
            if strcmp(sites.site_names(i), stable_names(j)) % if names are equal
                sites_occup(j) = sites_occup(j) + sites.occupancy(i);
                occup_parts(j,:) = occup_parts(j,:) + sites.occup_parts(i,:);
                site_count(j) = site_count(j) + 1;
                %fprintf('same: %s %s \n', sites.site_names{i}, stable_names{j})
            end
        end
    end
    
    % Fraction of their time atoms are at a given type of site:
    atom_loc_parts = zeros(nr_stable_names, sites.nr_parts);
    for j = 1:nr_stable_names
        atom_locations(j) = sites_occup(j)/(nr_steps*nr_atoms);
        for k = 1:sites.nr_parts
            atom_loc_parts(j,k) = occup_parts(j,k)/( ceil(nr_steps/sites.nr_parts) *nr_atoms);
        end
    end
    
    % Check how much time the atoms are not at 'known' positions
    total_frac = sum(atom_locations);
    fprintf('Fraction of time the diffusing atoms are at given sites: %d \n', ...
        total_frac)
    if total_frac < 0.75
        disp('WARNING! WARNING! WARNING!')
        disp('WARNING! Atoms spent more as 20% of the time at unknown locations!')
        disp('Change or add site locations for a better description of the jump and diffusion process!')
        disp('You can check if the sites are correct using the density plot (in plot_displacement.m)')
    end
    
    for j = 1:nr_stable_names
        sites_occup(j) = sites_occup(j)/(nr_steps*site_count(j));
        for k = 1:sites.nr_parts
            occup_parts(j,k) = occup_parts(j,k)/( ceil(nr_steps/sites.nr_parts) *site_count(j));
        end
    end

end