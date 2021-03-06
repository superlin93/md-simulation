function sim_data = read_vasp(outcar_file, vasprun_file, equil_time, diff_elem, output_file, diff_dim, z_ion)   
% Reads in an OUTCAR-file from VASP
% Made based on OUTCAR-files from VASP version 5.3.3 and version 5.4.1
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
    %% Read OUTCAR:
    sim_data = read_outcar(outcar_file, vasprun_file, sim_data);
    
    %% Calculate Attempt frequency, and standard deviation of vibration distance
    [sim_data.attempt_freq, sim_data.vibration_amp, sim_data.std_attempt_freq] = ...
        vibration_properties(sim_data, false);
    % Tracer diffusion coefficient and conductivity
    [sim_data.tracer_diffusion, sim_data.tracer_conductivity, sim_data.particle_density, ...
        sim_data.mol_per_liter] = tracer_properties(sim_data);
    % Save sim_data in a .mat file:
    save(output_file, 'sim_data'); 
    
    fprintf('Finished reading in the settings and atomic positions after %f minutes \n', toc/60 ) 
end

function sim_data = read_outcar(outcar_file, vasprun_file, sim_data)
    time = 0;
    nr_elements = 0;
    tic

%%  First read in things that are constant throughout the MD simulation:
    file = fopen(outcar_file);
    line = fgetl(file); %first line
    while time == 0
        temp = strsplit(line);
        if size(temp,2) > 1 % start comparing the text:
            if strcmp(temp{2}, 'POSITION') && strcmp(temp{3}, 'TOTAL-FORCE') % defines the start of the first time step
                time = time + 1; % Don't read in the positions at the first time step, because this is never necessary due to the equilibration time
            elseif strcmp(temp{2}, 'NSW') % The total number of MD steps
                total_steps = str2double(temp{4});	
            elseif strcmp(temp{2}, 'direct') % The crystal lattice
                for i = 1:3
                    line = fgetl(file); 
                    temp = strsplit(line);
                    for j = 2:4
                        sim_data.lattice(i,j-1) = str2double(temp{j});
                    end
                end
            elseif strcmp(temp{2}, 'energy') && strcmp(temp{4}, 'atom') %line before 'energy of atom' contains element names
                temp_prev = strsplit(line_prev);
                nr_elements = nr_elements + 1;
                sim_data.elements{nr_elements,1} = temp_prev{3}; % The element
            elseif strcmp(temp{2}, 'ions') % nr of ions per element
                for j = 6:5+nr_elements
                    per_elem = str2double(temp{j});
                    sim_data.nr_per_element(j-5, 1) = per_elem;
                end
                sim_data.nr_atoms = sum(sim_data.nr_per_element);
            elseif strcmp(temp{2}, 'TEBEG') % Temperature of the MD simulation
                temp = strsplit(line,{' ',';'}) ;% to remove the ;
                sim_data.temperature = str2double(temp{4});
            elseif strcmp(temp{2}, 'POTIM') % Size of the timestep (*1E-15 = in femtoseconds)
                sim_data.time_step = str2double(temp{4})*1E-15;
            elseif strcmp(temp{2}, 'volume') % Volume of the crystal lattice simulated
                sim_data.volume = str2double(temp{6})*1E-30; %Volume of the simulation, in m^3
            elseif strcmp(temp{2}, 'Mass') %Atomic masses as used by VASP (in u)
                line = fgetl(file);
                temp = strsplit(line);
                if length(temp) < 3+nr_elements %if mass is above 100 Vasp doesn't put a space between masses, so we have to split the string 
                    for j = 4:length(temp) %check where the error is
                        if length(temp{j}) > 6 %length of 6 is something like 100.99, so anything longer is wrong.
                            %extract the correct masses
                            correct_mass = strsplit(regexprep(temp{j},'\d*\.\d\d',' $0')); 
                            correct_mass = correct_mass(2:end); 
                            temp_mass = [temp, ' '];
                            if length(correct_mass) > 2
                                for k = 3:length(correct_mass); temp_mass = [temp_mass, ' ']; end;
                            end
                            temp_mass(j+2:end) = temp(j+1:end);
                            for k = 1:length(correct_mass)
                               temp_mass(j+k-1) = correct_mass(k); 
                            end
                            temp = temp_mass;
                        end
                    end
                end
                
                for j = 4:3+nr_elements
                    elem_mass = str2double(temp{j});
                    sim_data.element_mass(j-3) = elem_mass;
                end                
            end
        end
        line_prev = line; %the previous line
        line = fgetl(file);
    end  
    % Assign after reading in the constant things:
    sim_data.nr_elements = nr_elements;
    
    %Assign elements and diffusing element:
    counter = 1;
    sim_data.atom_element = cell(sim_data.nr_atoms,1);
    sim_data.diffusing_atom = false(sim_data.nr_atoms,1);
    
    sim_data.nr_diffusing = 0;
    for i = 1:sim_data.nr_elements
       % Where to find the diffusing atoms in sim_data, give the diffusing atoms the value 'true' 
        if strcmp(sim_data.elements{i}, sim_data.diff_elem)
            sim_data.nr_diffusing = sim_data.nr_diffusing + sim_data.nr_per_element(i);
            for j = counter:(counter + sim_data.nr_per_element(i) - 1) 
                sim_data.diffusing_atom(j) = true;
            end    
        end
        % Which element each atom is:
        for j = counter:(counter + sim_data.nr_per_element(i) - 1) 
            sim_data.atom_element{j} = sim_data.elements{i};
        end
        counter = counter + sim_data.nr_per_element(i);
    end
    
    % Check if the given diffusing element is present:
    if sim_data.nr_diffusing == 0
        error('ERROR! Given diffusing element not found in inputfile! Check your input')
    end
    
    %% Check and define simulation dependent things:
    sim_data.equilibration_steps = sim_data.equilibration_time/sim_data.time_step;
    % Number of steps to be used for the analysis:
    if isnan(total_steps)
        disp('WARNING! Total number of steps is undefined in OUTCAR, assuming 1 million steps!') %After values of 1 million VASP writes *******
        disp('WARNING! This will be adjusted after the atomic positions have been read in.')
        total_steps = 1000000; 
    end
    sim_data.nr_steps = round(total_steps - sim_data.equilibration_steps); 
    fprintf('Throwing away the first %4.0f steps because of the chosen equilibration time. \n', sim_data.equilibration_steps)
    
%% Now read positions of all atoms during the simulation:
    pos_line = ' POSITION                                       TOTAL-FORCE (eV/Angst)'; %KEEP THE SPACES!
    % Using the entire pos_line is approx. 4 times faster than splitting and checking with strcmp(temp{3}, 'TOTAL-FORCE')
    
    sim_data.cart_pos = zeros(3, sim_data.nr_atoms, sim_data.nr_steps); % Define cartesian positions array
    nr_atoms = sim_data.nr_atoms;
    step = 0;
    skip_steps = sim_data.equilibration_steps;
    line = fgetl(file);
    while ischar(line)   
        if strcmp(pos_line, line)
            time = time + 1;
            if time > skip_steps % ! Equilibration steps are thrown away!
                step = step + 1; % The time step 
                fgetl(file); %skip a line
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

    if step == 0 || time < 0.25*sim_data.nr_steps % No or very few positions found in OUTCAR, so read them from vasprun.xml
        disp('WARNING! The OUTCAR-file is missing a lot of atomic positions from the MD simulation!')
        if ~exist(vasprun_file, 'file')
            error('ERROR! Put vasprun.xml in the folder to read the atomic positions from that file, then run analyse_md again')
        else % USE vasprun.xml to read in coordinates
            [sim_data, step] = read_vasprunxml(vasprun_file, sim_data);
        end
    end
    
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

function [sim_data, step] = read_vasprunxml(vasprun_file, sim_data)
% Reads the atomic positions during the MD simulation from vasprun.xml 
% only when the atomic coordinates are missing in OUTCAR! (for Shiv)
    disp('WARNING! The atomic positions during the MD simulation are read from vasprun.xml')
    %Start:
    file = fopen(vasprun_file);
    pos_line = '   <varray name="positions" >';
    nr_atoms = sim_data.nr_atoms;
    step = 0;
    lattice = sim_data.lattice;
    skip_steps = sim_data.equilibration_steps;
    frac_pos = zeros(3,1);
    time = 0;
    
    tic
    line = fgetl(file); %the first line
    while ischar(line)
        if strcmp(pos_line, line)
            time = time + 1;
            if time > skip_steps
                step = step + 1;
                %Faster way should be possible, by reading all coor at
                %once instead of the loops, but good enough for now
                for atom = 1:nr_atoms %loop over the atoms              
                    line = strsplit(fgetl(file)); %the next line
                    for j = 1:3 
                        frac_pos(j) = sscanf(line{j+2}, '%f');
                    end
                    % Not very efficient, but easy to implement:
                    cart = frac_to_cart(frac_pos, lattice);
                    sim_data.cart_pos(:,atom,step) = cart;
                end
            end % < time
            if mod(step+1,2500) == 0 %To see that stuff is still happening:
                fprintf('Reading timestep %d after %f seconds. \n', step, toc)
            end
        end %strcmp
        line = fgetl(file); %the next line         
    end
    fclose(file);
end