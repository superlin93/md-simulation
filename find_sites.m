function [sites, finished] = find_sites(sim_data, material, nr_parts)
% Given the locations of sites determine at which site the atoms are.
% Materials known: 'argyrodite', 'latp', 'LiSnPS', 'na3ps4' and
% 'li3ps4_beta', 'MnO2_lambda'
    disp('Finding site occupations for given structure')
    finished = false; % Check if the function has finished correctly
%% Get the names and positions for the given material:
    [names, sites_pos, supercell] = known_materials(material);
    sites.names = unique(names);
    sites.supercell = supercell;
    
    %% Define some things:
    dist_close = 2*sim_data.vibration_amp; % 'amplitude of vibration'
    nr_sites = size(sites_pos,2); 
    part_size = ceil(sim_data.nr_steps/nr_parts);   
    lattice = sim_data.lattice;

    %% The cartesian positions of the given sites:
    cart_pos = zeros(3, nr_sites);
    for i = 1:nr_sites
        cart_pos(:,i) = frac_to_cart(sites_pos(:,i), lattice);
    end
       
    %% Names of the sites:
    nr_per_name = zeros(size(sites.names, 1),1);
    for i = 1:size(names, 1)
        for j = 1:size(sites.names, 1)
            if strcmp(sites.names{j}, names{i})
               nr_per_name(j) = nr_per_name(j) + 1;
            end
        end
    end
    sites.nr_per_name = nr_per_name;
    
 %% Add the names of transition sites: 
    unique_names = unique(names);
    % Single names:
    for i = 1:size(unique(names))
        trans_name = strcat('Transition:', unique_names(i));
        a = length(sites.names) + 1;
        sites.names(a) = trans_name;
    end
    % All combinations:
    for i = 1:size(unique(names))
        for j = 1:size(unique(names))
            trans_name = strcat('Transition:', unique_names(i), '-', unique_names(j));
            a = length(sites.names) + 1;
            sites.names(a) = trans_name;
        end
    end

%% Loop over all possible positions, to see how close they are to each other
% Should be less as the 2*dist_close to prevent overlap, otherwise it will be adjusted
    for site1 = 1:nr_sites-1
        pos1 = sites_pos(:,site1);
        for site2 = site1+1:nr_sites
            pos2 = sites_pos(:,site2);
            dist = sqrt(calc_dist_sqrd_frac(pos1, pos2, lattice));
            if dist < 2*dist_close % Make dist_close smaller to prevent overlapping sites:
                fprintf('WARNING! Crystallographic sites are overlapping with the chosen dist_close (%d), making dist_close smaller! \n', dist_close)
                dist_close = 0.5*(dist - 0.01); % -0.01 for extra certainty
                fprintf('dist_close changed to %d \n', dist_close)
                if dist_close < 0.25 %Angstrom
                    disp('ERROR! ERROR! Two crystallographic sites are within half an Angstrom of each other')
                    disp('ERROR! This is NOT realistic, check/change the given crystallographic site locations in known_materials.m!')
                    fprintf('ERROR! Sites number %d and %d are overlapping \n', site1, site2) 
                    fprintf('ERROR! Coordinates are: (%d, %d, %d) and (%d, %d, %d) \n', ... 
                    pos1(1), pos1(2), pos1(3), pos2(1), pos2(2), pos2(3))                    
                    return
                end
            end
        end
    end
    close_sqrd = dist_close*dist_close; % to avoid taking the root everytime
% Necessary matrices:
    atom_site = zeros(sim_data.nr_steps, sim_data.nr_diffusing); % A value of zero means its not at a known site
    site_occup = zeros(nr_sites,nr_parts);
    transitions = zeros(nr_sites);
    trans_counter = 0;
    all_trans = zeros(1,4);
    succes = zeros(nr_sites, nr_sites, nr_parts);
%%  Loop over all times and diffusing atoms: 
    diff_atom = 0;
    part = 1;
    end_step = sim_data.nr_steps;
    disp('Determining the site positions of all diffusing atoms:')
    for atom = sim_data.start_diff_elem:sim_data.end_diff_elem
        diff_atom = diff_atom + 1;
        % First find the first site this atom at:
        not_found = true;
        site = 0;
        time_start = 1;
        while not_found && time_start < end_step
            if site < nr_sites
                site = site + 1;
            else
                time_start = time_start + 1;
                site = 1;
            end
            pos_atom = sim_data.frac_pos(:,atom,time_start);             
            dist = calc_dist_sqrd_frac(pos_atom, sites_pos(:, site), lattice);
            if dist < close_sqrd
                not_found = false;
                % Atom is at this site
                atom_site(time_start, diff_atom) = site;
                % Increase occupancy of this site    
                site_occup(site, part) = site_occup(site, part) + 1; 
                prev_site = site; % Update previous site for next transition
            end
        end               
        % After the first site has been found:
        for time = time_start+1:sim_data.nr_steps   
           % Divide the simulation in multiple part's for statistics:
            part = ceil(time/part_size);
            pos_atom = sim_data.frac_pos(:,atom,time);           
        % First check the site found at the last time_step, since it's the most likely position:
            dist = calc_dist_sqrd_frac(pos_atom, sites_pos(:, prev_site), lattice);
            if dist < close_sqrd
                atom_site(time, diff_atom) = prev_site; % Atom is at this site   
                site_occup(prev_site, part) = site_occup(prev_site, part) + 1; % Increase occupancy of this site
                trans_start = time + 1;
            else %not found at the previous position: 
            % Loop over all sites, exit the loop once a site has been found.
                not_found = true;
                site = 0;
                while not_found && site < nr_sites
                    site = site + 1;
                    dist = calc_dist_sqrd_frac(pos_atom, sites_pos(:, site), lattice);
                    if dist < close_sqrd
                        not_found = false;
                    % Atom is at this site
                        atom_site(time, diff_atom) = site;
                    % Increase occupancy of this site    
                        site_occup(site, part) = site_occup(site, part) + 1;    
                    end
                end
            % If a transition happened, remember some things:
                if ~not_found %Transition happened, i.e. the atom is found at a new site
                    transitions(prev_site, site) = transitions(prev_site, site) + 1;
                    trans_counter = trans_counter + 1;
                    all_trans(trans_counter,1) = diff_atom; %the atom
                    all_trans(trans_counter,2) = prev_site; % the starting position
                    all_trans(trans_counter,3) = site; %the final position     
                    all_trans(trans_counter,4) = trans_start; %the time the atom left it's previouw position
                    all_trans(trans_counter,5) = time; % the timestep (at the end of the transition)
                    succes(prev_site, site, part) = succes(prev_site, site, part) + 1;  
                    prev_site = site;
                end
            end
        end  
        fprintf('*') % One atom has been finished...
    end
    fprintf('\n') % go to next line in output
    fprintf('Number of transitions between sites: %d \n', trans_counter)
    
    % Also assign the (failed) transition states with a site-number and name
    diff_atom = 0;
    time_start = 0;
    nr_trans_found = 0;
    for atom = sim_data.start_diff_elem:sim_data.end_diff_elem
        diff_atom = diff_atom + 1;
        prev_site = 0;
        % What to do with atoms starting outside a given site?
        for time = time_start+1:sim_data.nr_steps
            if atom_site(time, diff_atom) == 0 && prev_site == 0
                % It left the site
                if time ~= 1 
                    prev_site = atom_site(time - 1, diff_atom);
                    trans_start = time;
                end
            elseif atom_site(time,diff_atom) ~= 0 && prev_site ~= 0   
                % It appeared at a new site
                now_site = atom_site(time, diff_atom);
                trans_end = time-1;
                % Determine the 'trans_site_nr'
                if now_site == prev_site
                    trans_site_nr = nr_sites + now_site;
                    names(trans_site_nr) = strcat('Transition:', names(now_site));
                    atom_site(trans_start:trans_end, diff_atom) = trans_site_nr;
                    prev_site = 0;
                else % loop over all transitions already found
                    if nr_trans_found == 0 % For the first transition site:
                       trans_found(1,1) = prev_site;
                       trans_found(1,2) = now_site;
                       nr_trans_found = 1;
                       trans_nr = 1;
                    else
                        trans_nr = 0;
                        % Check if it has been found earlier:
                        for i = 1:nr_trans_found
                            if trans_found(i,1) == prev_site && trans_found(i,2) == now_site
                                trans_nr = i;
                            end
                        end
                        if trans_nr == 0 % Not found in previous step, so it is a NEW one
                            nr_trans_found = nr_trans_found + 1;
                            trans_found(nr_trans_found,1) = prev_site;
                            trans_found(nr_trans_found,2) = now_site;
                            trans_nr = nr_trans_found;
                        end
                    end
                    % Assign trans_site_nr etc.:
                    trans_site_nr = 2*nr_sites + trans_nr;
                    atom_site(trans_start:trans_end, diff_atom) = trans_site_nr;
                    names(trans_site_nr) = strcat('Transition:', names(prev_site), '-', names(now_site));
                    prev_site = 0;
                end
            end
        end
    end
    % Save things to sites:
    sites.site_radius = dist_close;
    sites.frac_pos = sites_pos;
    sites.cart_pos = cart_pos;
    sites.site_names = names;
    sites.occupancy = sum(site_occup, 2);
    sites.occup_parts = site_occup;
    sites.atoms = atom_site;
    sites.transitions = transitions;
    sites.all_trans = all_trans;
    sites.succes = succes;
    sites.nr_parts = nr_parts;
    
    finished = true;
end    