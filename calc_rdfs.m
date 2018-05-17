function rdf = calc_rdfs(sites, sim_data, res, max_dist)
% Calculate the radial distribution function for the diffusing atom (overall and per site)
    disp('Calculating Radial Distribution Functions for the diffusing element.... ')
    %tic
    max_bin = ceil(max_dist/res);
    fprintf('Maximum distance for RDF: %f Angstrom \n', max_dist)
    fprintf('Resolution for RDF: %f Angstrom \n', res)   
    
    nr_atoms = sim_data.nr_atoms;
    nr_elem = sim_data.nr_elements;
    lat = sim_data.lattice;

    nr_rdfs = numel(sites.names); %nr of different names
    rdf_names = sites.names;

    % Use elem_struct for faster searching:
    elem_struct = zeros(nr_elem+1,1);
    elem_struct(1) = 1;
    for i = 2:nr_elem+1
        elem_struct(i) = elem_struct(i-1) + sim_data.nr_per_element(i-1);
    end
    
    % Sizes of all the necessary things:
    rdfs = zeros(nr_elem, max_bin, nr_rdfs);
    integrated = zeros(nr_elem, max_bin, nr_rdfs);
    norm_counter = zeros(nr_rdfs,1);
    all = zeros(nr_elem, max_bin);
    atom_rdf = zeros(nr_elem, max_bin);
%% Main loop:    
    diff_elem = 0;
    for atom = sim_data.start_diff_elem:sim_data.end_diff_elem %Only for the diffusing element
        diff_elem = diff_elem + 1;
        for time = 1:sim_data.nr_steps
            atom_rdf(:,:) = 0.0; % Calculate RDF for THIS atom at THIS time step
            pos1 = sim_data.frac_pos(:, atom, time); %atom position
            for i = 1:nr_atoms % Calculate distance from 'atom' to all other atoms
                if i ~= atom % But not with itself
                    dist = sqrt(calc_dist_sqrd_frac(pos1, sim_data.frac_pos(:, i, time), lat));
                    dist_bin = ceil(dist/res);
                    if dist_bin <= max_bin 
                        % Determine which element the 'other atom' is
                        elem = 0;
                        j = 0;
                        while elem == 0
                            j = j + 1;
                            % Search for element:
                            if i >= elem_struct(j) && i < elem_struct(j+1)
                                elem = j;
                            end
                        end
                        atom_rdf(elem, dist_bin) = atom_rdf(elem, dist_bin) + 1;
                    end
                end
            end %i
            % RDF determined for this atom at this time:
            all = all + atom_rdf; % all = RDF for the entire simulation
            if sites.atoms(time, diff_elem) ~= 0 % FOR NOW A QUICK FIX, should never be zero
                site_name = sites.site_names(sites.atoms(time, diff_elem)); %sites.atoms only contain diff_elem info
                name = 0;
                not_found = true;
                while not_found % Find the nr. corresponding to the correct site name
                    name = name +1;
                    if strcmp(site_name, rdf_names{name})
                        rdfs(:,:,name) = rdfs(:,:,name) + atom_rdf; % Add RDF to the correct name
                        norm_counter(name) = norm_counter(name) + 1; % Normalisation factor = the nr. of timesteps it is occupied
                        not_found = false;
                    end
                end
            end 
        end
        fprintf('*') % Shows that one atom is finished
    end
    fprintf(' \n')
    %fprintf('Finished calculating RDFs after %f minutes. \n', toc/60)
    
    % Integrate the RDF's
    for i = 1:nr_rdfs
        for j = 2:max_bin % Should be ok since the first bin has a distance of 0.5*res, so should always be empty.
            integrated(:,j,i) = integrated(:,j-1,i) + rdfs(:,j,i);
        end
    end
    
    % create data structure:
    rdf.distributions = rdfs;
    rdf.integrated = integrated;
    rdf.rdf_names = rdf_names;
    rdf.elements = sim_data.elements;
    rdf.max_dist = max_dist;
    rdf.resolution = res;
    rdf.total = all;
end