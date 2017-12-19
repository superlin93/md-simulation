function sim_data = read_vasp(outcar_file, equil_time, diff_elem, output_file, diff_dim, z_ion)   
% Reads in an OUTCAR-file from VASP
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
    sim_data = read_outcar(outcar_file, sim_data);
    
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

function sim_data = read_outcar(outcar_file, sim_data)
    time = 0;
    tic
    pos_line = ' POSITION                                       TOTAL-FORCE (eV/Angst)'; %KEEP THE SPACES!
%%  First read in things that are constant throughout the MD simulation:
    file = fopen(outcar_file);
    line = fgetl(file); %first line
    while time == 0
        temp = strsplit(line);
        if size(temp,2) > 1 % start comparing the text:
            if strcmp(temp{2}, 'POSITION') % defines the start of a new time step
                time = time + 1;
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
            elseif strcmp(temp{2}, 'INCAR:') % the number of elements
                nr_elements = 0;
                line = fgetl(file);
                temp = strsplit(line);
                while size(temp,2) > 1 && strcmp(temp{2}, 'POTCAR:') 
                    nr_elements = nr_elements + 1;
                    sim_data.elements{nr_elements,1} = temp{4}; % The element
                    line = fgetl(file);
                    temp = strsplit(line);          
                end
                sim_data.nr_elements = nr_elements;
            elseif strcmp(temp{2}, 'ions') % nr of ions per element
                for j = 6:5+sim_data.nr_elements
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
                temp = strsplit(fgetl(file));
                for j = 4:3+sim_data.nr_elements
                    elem_mass = str2double(temp{j});
                    sim_data.element_mass(j-3) = elem_mass;
                end                
            end
        end
        line = fgetl(file);
    end
    %% Check and define simulation dependent things:
    sim_data.equilibration_steps = sim_data.equilibration_time/sim_data.time_step;
    % Number of steps to be used for the analysis:
    sim_data.nr_steps = round(total_steps - sim_data.equilibration_steps); 
    fprintf('Throwing away the first %4.0f steps because of the chosen equilibration time. \n', sim_data.equilibration_steps)

    % Define cartesian positions array
    sim_data.cart_pos = zeros(3, sim_data.nr_atoms, sim_data.nr_steps); 

    % Diffusing element specific:
    counter = 1;
    sim_data.atom_element = cell(sim_data.nr_atoms,1);
    sim_data.nr_diffusing = 0;
    for i = 1:sim_data.nr_elements
        if strcmp(sim_data.elements{i}, sim_data.diff_elem)
            sim_data.nr_diffusing = sim_data.nr_per_element(i);
            % Where to find the diffusing atoms in sim_data.cart_pos (and others):
            sim_data.start_diff_elem = counter; 
            sim_data.end_diff_elem = counter + sim_data.nr_per_element(i) -1;
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
    
%% Now read positions of all atoms during the simulation:
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
                fprintf('Reading timestep %d after %f minutes. \n', step+1, toc/60)
            end
        end %positions
        line = fgetl(file); %the next line         
    end
    fclose(file);

%% Reading OUTCAR file done:
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