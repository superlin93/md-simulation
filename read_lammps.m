function sim_data = read_lammps(xyz_file, lammps_input, lammps_output, lammps_structure, ...
                 equil_time, diff_elem, output_file, diff_dim, z_ion)
% Reads in data from LAMMPS simulations
    %% Define constants
    sim_data.e_charge = 1.60217657e-19; %Electron charge
    sim_data.k_boltzmann = 1.3806488e-23; %Boltzmann's constant 
    sim_data.avogadro = 6.022140857e23; %Avogadro's number
 
    sim_data.diffusion_dim = diff_dim;
    sim_data.ion_charge = z_ion;
    
    %% Initialise:    
    sim_data.lattice = zeros(3);
    sim_data.diff_elem = diff_elem;
    sim_data.equilibration_time = equil_time;
    
    %% Read output of  simulation:
    sim_data = read_lammpsfiles(xyz_file, lammps_input, lammps_output, lammps_structure, sim_data);
    
    %% Calculate Attempt frequency, and standard deviation of vibration distance
    [sim_data.attempt_freq, sim_data.vibration_amp, sim_data.std_attempt_freq] = ...
        vibration_properties(sim_data, false);
    % Tracer diffusion coefficient and conductivity
    [sim_data.tracer_diffusion, sim_data.tracer_conductivity, sim_data.particle_density, ...
        sim_data.mol_per_liter] = tracer_properties(sim_data);
    % Save sim_data in a .mat file:
    save(output_file, 'sim_data'); 
    
    fprintf('Finished reading in the OUTCAR file after %f minutes \n', toc/60 ) 
end

function sim_data = read_lammpsfiles(xyz_file, lammps_input, lammps_output, lammps_structure, sim_data)
     tic
%%  First read in things that are constant throughout the MD simulation:
    % From the input file of LAMMPS get some settings:
     file = fopen(lammps_input);
     line = fgetl(file); %first line
     while ischar(line)
        temp = strsplit(line);
        if size(temp,2) > 1 % start comparing the text:     
            if strcmp(temp{1}, 'run') % The total number of MD steps    
                total_steps = str2double(temp{2});
            elseif strcmp(temp{2}, 'T1') && strcmp(temp{3}, 'equal') % Temperature of the MD simulation
                sim_data.temperature = str2double(temp{4});
            elseif strcmp(temp{1}, 'timestep') % Size of the timestep (*1E-12 = in picoseconds)
                sim_data.time_step = str2double(temp{2})*1E-12;   
            end
        end
        line = fgetl(file);
     end
     fclose(file);
     
     % From the outputfile of LAMMPS get the (optimised) volume and lattice
     file = fopen(lammps_output);
     line = fgetl(file); %first line
     no_tilt = true;
     while ischar(line)
        temp = strsplit(line);
        if size(temp,2) > 1 % start comparing the text: 
         % The crystal lattice
            if strcmp(temp{2}, 'v_Timer')
                line_before = temp; % remember what's on the line before
                
                line = fgetl(file); % next line contains the numbers
                temp = strsplit(line);
                sim_data.lattice(1,1) = str2double(temp{8});
                sim_data.lattice(2,2) = str2double(temp{9});
                sim_data.lattice(3,3) = str2double(temp{10});
                
                if strcmp(line_before{10}, 'Xy') % For non-rectangular lattices
                    no_tilt = false;
                    sim_data.lattice(2,1) = str2double(temp{11});
                    sim_data.lattice(3,1) = str2double(temp{12});
                    sim_data.lattice(3,2) = str2double(temp{13});
                    sim_data.volume = str2double(temp{14})*1E-30; %Volume of the simulation, in m^3    
                else % Assume a rectangular lattice as there is no tilt given...
                    sim_data.volume = str2double(temp{11})*1E-30; %Volume of the simulation, in m^3  
                end
             end
         end
         line = fgetl(file);
     end
     
     if no_tilt
        disp('WARNING! No tilt found in the LAMMPS output')
        disp('WARNING! Assuming a rectangular lattice!')
     end
     fclose(file);
     
    % From the structure file of LAMMPS:
    file = fopen(lammps_structure);
    line = fgetl(file); %first line
     
     while ischar(line)
        temp = strsplit(line);
        if size(temp,2) > 1
            if strcmp(temp{2}, 'atoms') 
                sim_data.nr_atoms = str2double(temp{1}); %total number of atoms in the simulation
                sim_data.atom_element = cell(sim_data.nr_atoms,1);
            elseif size(temp,2) > 2 %Amount of different elements
                if strcmp(temp{2}, 'atom') && strcmp(temp{3}, 'types')
                    sim_data.nr_elements = str2double(temp{1});
                    % Define for later:
                    sim_data.nr_per_element = zeros(sim_data.nr_elements,1);
                    sim_data.elements = cell(sim_data.nr_elements,1);
                end
            end
        elseif strcmp(temp{1}, 'Masses')
            line = fgetl(file); %skip this line
            for i = 1:sim_data.nr_elements
                line = fgetl(file);
                temp = strsplit(line);
                sim_data.element_mass(i) = str2double(temp{2});
                % Get the element name based on the mass:
                sim_data.elements{i} = element_from_mass(sim_data.element_mass(i));
            end
        elseif strcmp(temp{1}, 'Atoms') %Start counting the number of atoms per element
            line = fgetl(file); %skip this line
            for j = 1:sim_data.nr_atoms 
                line = fgetl(file);
                temp = strsplit(line);
                if strcmp(temp{1}, '') %Sometimes temp starts with '', sometimes it does not...
                    nr = str2double(temp{4});
                else 
                    nr = str2double(temp{3}); 
                end
                % count the number of atoms per element:
                sim_data.nr_per_element(nr) = sim_data.nr_per_element(nr) + 1;
                % remember which element each atom is:
                sim_data.atom_element{j} = sim_data.elements{nr};
            end
        end    
        line = fgetl(file);
     end
     fclose(file);
     
     for i = 1:sim_data.nr_elements
         fprintf('Found %4.0f atoms of element %3s \n', sim_data.nr_per_element(i), sim_data.elements{i})
     end          
     
    % Diffusing element specific:
    counter = 1;
    sim_data.nr_diffusing = 0;
    for i = 1:sim_data.nr_elements
        if strcmp(sim_data.elements{i}, sim_data.diff_elem)
            sim_data.nr_diffusing = sim_data.nr_per_element(i);
            % Where to find the diffusing atoms in sim_data.cart_pos (and others):
            sim_data.start_diff_elem = counter; 
            sim_data.end_diff_elem = counter + sim_data.nr_per_element(i) -1;
        end
    end
    
    % Check if the given diffusing element is present:
    if sim_data.nr_diffusing == 0
        error('ERROR! Given diffusing element not found in inputfile! Check your input')
    end
    
    %% Check and define simulation dependent things:
    sim_data.equilibration_steps = sim_data.equilibration_time/sim_data.time_step;
    % Number of steps to be used for the analysis:
    sim_data.nr_steps = round(total_steps - sim_data.equilibration_steps); 
    fprintf('Throwing away the first %4.0f steps because of the chosen equilibration time. \n', sim_data.equilibration_steps)

%% Now read positions of all atoms during the simulation from the xyz-file:
    file = fopen(xyz_file);

    sim_data.cart_pos = zeros(3, sim_data.nr_atoms, sim_data.nr_steps); % Define cartesian positions array
% For testing avoid arrays of several gigabytes:
%    sim_data.cart_pos = zeros(3, sim_data.nr_atoms, 50);
    nr_atoms = sim_data.nr_atoms;
    
    pos_line = 'Atoms.'; %after which word the positions of the next timestep begin
    
    time = 0;
    step = 0;
    skip_steps = sim_data.equilibration_steps;
    line = fgetl(file);
    while ischar(line)
        check_start = strsplit(line);
        if strcmp(pos_line, check_start(1) ) % Check if the line starts with pos_line           
            time = time + 1;
            if time > skip_steps % ! Equilibration steps are thrown away!
                step = step + 1; % The time step 
                for atom = 1:nr_atoms %loop over the atoms     
                    line = strsplit(fgetl(file)); %the next line
                    for j = 1:3 
                        sim_data.cart_pos(j,atom,step) = sscanf(line{j+1}, '%f');
                    end
                end
            end
            if mod(step+1,2500) == 0 % Show that reading in is still happening:
                fprintf('Reading timestep %d of %d after %f minutes. \n', ... 
                    step+1, sim_data.nr_steps, toc/60)
            end
        end %positions
        line = fgetl(file); %the next line         
    end
    fclose(file);
    
%% Reading positions done:
    % If the simulation was not completely finished:
    if step ~= sim_data.nr_steps
        sim_data.nr_steps = step;
        temp_cart = sim_data.cart_pos(:, :, 1:step);
        sim_data.cart_pos = zeros(3, sim_data.nr_atoms, sim_data.nr_steps);
        sim_data.cart_pos = temp_cart;
    end
    % Total simulated time:
    sim_data.total_time = sim_data.nr_steps*sim_data.time_step;
    
    % Determine fractional positions and displacement:
    sim_data = frac_and_disp(sim_data); 
    
end

function element = element_from_mass(mass)
% Determine the name of the element based on the given mass:
    if mass == 1
        element = 'H';
    elseif mass == 7
        element = 'Li';
    elseif mass == 16
        element = 'O';
    elseif mass == 23
        element = 'Na';        
    elseif mass == 27
        element = 'Al';        
    elseif mass == 31
        element = 'P';
    elseif mass == 32
        element = 'S';        
    else
        error('Element not recognised based on its mass! Add it to element_from_mass in read_lammps.m')
    end
end