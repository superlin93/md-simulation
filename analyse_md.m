function analyse_md(folder, diff_elem, material)
% Materials known: 'argyrodite', 'latp', 'LiSnPS', 'na3ps4', 'li3ps4_beta' and 'MnO2_lambda'

% When you use this code in academic work please cite the accompanying paper: 
% Analysis of diffusion in solid state electrolytes through MD-simulations, 
% improvement of the Li-ion conductivity in \ce{\beta-Li3PS4} as an example;
% Niek J.J. de Klerk, Eveline van der Maas and Marnix Wagemaker (submitted)

%% Settings:
% !!! WARNING! The settings for equil_time, diffusion_dimensions, z_ion,
% !!! nr_parts and dist_collective are ONLY used when creating the sim_data.mat and sites.mat files !!!
% !!! RENAME OR REMOVE THESE FILES IF YOU CHANGE THESE SETTINGS !!!   
    equil_time = 2.5E-12; % Equilibration time in seconds (1E-12 = 1 picosecond)  
    diffusion_dimensions = 3; % number of diffusion dimensions
    z_ion = 1.0; % Ionic charge of the diffusing ion
    nr_parts = 10; % In how many parts to divide your simulation for statistics:   
    dist_collective = 4.5; % Maximum distance for collective motions in Angstrom
    
    % Pictures 
    show_pics = true; % Show plots (or not)
    density_resolution = 0.2; % Resolution for the diffusing atom density plot, in Angstrom  
    jump_res = 0.1; %Resolution for the nr. of jumps vs. distance plot,in Angstrom
    
    % Radial Distribution Functions
    % !!! WARNING! rdf_res and rdf_max_dist are ONLY used when creating the rdfs.mat file !!!    
    % !!! RENAME OR REMOVE THE RDF-FILE IF YOU WANT TO CHANGE THESE SETTINGS !!!   
    rdfs = true; % Calculate and show Radial Distribution Functions
    rdf_res = 0.1; % Resolution of the RDF bins in Angstrom
    rdf_max_dist = 10; %Maximal distance of the RDF in Angstrom 

    % Movie showing the jumps:  
    movie = false; % Make a movie showing the jumps (or not) 
    nr_steps_frame = 5; % How many time steps per frame in the movie, increase to get smaller and quicker movies
    start_end = [5000; 7500]; % The time-steps for which to make a movie ([start at step; end at step])
       
 %% Standard file names:
    outcar_file = [folder, '/OUTCAR'];
    vasprun_file = [folder, '/vasprun.xml']; %Backup for if the atomic positions are not in the OUTCAR
    sim_data_file = [folder, '/simulation_data.mat'];
    sites_file = [folder, '/sites.mat'];
    rdf_file = [folder, '/rdfs.mat'];
    movie_file = [folder, '/jumps_movie.mp4'];

%% Get simulation data:
    fprintf('Folder given is: %s \n', folder) 
    if ~exist(folder, 'dir') % Check if given folder exists
            disp('ERROR! Given folder does not exist!')
            return
    end
            
    % Check if simulation_data.mat exists in the given folder:    
    if ~exist(sim_data_file, 'file')
    	fprintf('%s not found. \n', sim_data_file) 
        % Check for OUTCAR file:
        if exist(outcar_file, 'file')
            % Read in the simulation data from outcar_file
            fprintf('Reading simulation data from %s, this can take a while.... \n', outcar_file)            
            sim_data = read_vasp(outcar_file, vasprun_file, equil_time, diff_elem, sim_data_file, diffusion_dimensions, z_ion);
        else
            fprintf('ERROR! %s not found! \n', outcar_file)   
            return
        end
    else % sim_data exists already:
        disp('Found simulation data file in given folder')
        load(sim_data_file) 
    end    
    
%% Find sites and transitions  
    % Check if sites_file already exists:
    if ~exist(sites_file, 'file')
        % Find the sites if the sites_file does not exist:
        [sites, finished]= find_sites(sim_data, material, nr_parts);
        if ~finished
            disp('Find sites exited with an error, not analysing this simulation further... ')
            return            
        end
    % Save the material
        sites.material = material;
    % Determine the fractional occupancies of sites:    
        [sites.stable_names, sites.sites_occup, sites.atom_locations, ... 
            sites.occup_parts, sites.atom_loc_parts] = calc_site_occups(sites);       
    % Calculate jump rates and activation energy for all the differently named 
    % transitions between the given sites:  
        [sites.jump_names, sites.nr_jumps, sites.rates, sites.e_act] = calc_rates(sites, sim_data);
    % Possible correlations between jumps
        [sites.collective, sites.coll_jumps, sites.coll_matrix, ...
            sites.multi_coll, sites.solo_jumps] = possible_collective(sites, sim_data, dist_collective);         
    % The fraction of solo jumps:  
        sites.solo_frac = sites.solo_jumps/sum(sites.nr_jumps);
    % The jump diffusivity:
        sites.jump_diffusivity = jumps_vs_dist(sites, sim_data, false, jump_res);
    % The correlation factor:
        sites.correlation_factor = sim_data.tracer_diffusion/sites.jump_diffusivity;
        save(sites_file, 'sites')
    else
        load(sites_file)
        fprintf('Found sites file: %s .\n', sites_file)
        if ~strcmp(sites.material, material)
            fprintf('ERROR! The material in sites.mat (%s) differs from the given one (%s)! \n', sites.material, material)
            fprintf('If this is not a mistake rename or remove sites.mat and try again. \n')
            return 
        end
    end
    
%% Show plots:
    if show_pics 
        % Plots based on atomic displacement/position:
        % Displacement per element, for diffusing element per atom, histogram, and density plot:
        plots_from_displacement(sim_data, sites, density_resolution); 
        % Jump sites:  
        plot_sites(sites, sim_data.lattice, true); %true is to show names of sites instead of numbers
        % Nr. of jumps vs. jump distance:
        jumps_vs_dist(sites, sim_data, true, jump_res);
        % Attempt frequency and vibration amplitude:
        vibration_properties(sim_data, true); % true is to show the figures
        % Collective motions:
        if sites.nr_jumps > 1
            plot_collective(sites);
        end
    end
    
    %% Radial distribution functions
    % Check if rdfs == true, and if the file does not exist yet:
    if rdfs && ~exist(rdf_file, 'file') 
        % rdf's per site
        fprintf('Radial Distribution Functions file %s not found \n', rdf_file)
        rdf = calc_rdfs(sites, sim_data, rdf_res, rdf_max_dist);
        save(rdf_file, 'rdf')
        plot_rdfs(rdf)
    elseif rdfs
        fprintf('Using Radial Distribution Functions file %s \n', rdf_file)
        load(rdf_file)
        plot_rdfs(rdf)
    end

    %% Movie of the transitions: 
    if movie && sum(sites.nr_jumps) > 0 %0 jumps will make a very boring movie
        if ~exist(movie_file, 'file')
            make_movie(sites, sim_data, start_end, movie_file, nr_steps_frame); 
        elseif exist(movie_file, 'file')
            fprintf('Movie file %s found, rename the movie file if you want to make another movie. \n', movie_file)
        end
    elseif movie && sum(sites.nr_jumps) == 0 %0 jumps will make a very boring movie
        fprintf('Not making a movie because there are no jumps occuring in the MD simulation. \n')
    end
    
    disp('Analysis of MD simulation done')
end
