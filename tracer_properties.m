function [tracer_diff, tracer_conduc, particle_density, mol_per_liter] = ... 
    tracer_properties(sim_data)
% Calculate the diffusion and conductivity based on the displacement
% For calculating the diffusion:
    ang2m = 1E-10; %converts Angstroms to meters
    dimensions = sim_data.diffusion_dim;
    z_ion = sim_data.ion_charge;
    
    disp('-------------------------------------------')
    disp('Calculating tracer diffusion and conductivity based on:')
    fprintf('%f dimensional diffusion, and an ion with a charge of %f \n', dimensions, z_ion)

    % The particle density (in particles/m^3)    
    particle_density = sim_data.nr_diffusing/sim_data.volume; %in particles/(m^3)
    % Concentration of diffusing atom in mol/L (1000 converts from m^3 to liters)
    mol_per_liter = particle_density/(1000*sim_data.avogadro);
    
    % Mean Squared Displacement:        
    displacement = sim_data.displacement(sim_data.start_diff_elem:sim_data.end_diff_elem, sim_data.nr_steps);
    sqrd_disp = displacement.^2;
    msd = mean(sqrd_disp); % In Angstrom^2

    % Diffusivity = MSD/(2*dimensions*time)
    tracer_diff = (msd * ang2m^2)/(2*dimensions*sim_data.total_time); % In m^2/sec

    % Conductivity = elementary_charge^2 * charge_ion^2 * diffusivity * particle_density / (k_B * T)
    tracer_conduc = ((sim_data.e_charge^2) * (z_ion^2) * tracer_diff * particle_density)/ ...
        (sim_data.k_boltzmann*sim_data.temperature); % In Siemens/meter

    fprintf('Tracer diffusivity determined to be (in meter^2/sec): %d \n', tracer_diff)
    fprintf('Tracer conductivity determined to be (in Siemens/meter): %d \n', tracer_conduc)
    disp('-------------------------------------------')
end