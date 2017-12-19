function jump_diff = jumps_vs_dist(sites, sim_data, pics, resolution)
% Calculate the jump diffusivity and plot the nr. of jumps vs. distance of the jumps
    nr_sites = size(sites.occupancy,1);
    nr_diffusing = size(sites.atoms,2);    
    ang2m2 = 1E-20; %Angstrom^2 to meter^2
    
    lattice = sim_data.lattice;
    dimensions = sim_data.diffusion_dim;
    sim_time = sim_data.total_time;
    
    max_dist = 10.0; %Much larger as the largest jump distance in typical materials
    jumps_vs_dist = zeros(ceil(max_dist/resolution),1);
    max_bin = 0;
    jump_diff = 0;
    for i = 1:nr_sites-1
        for j = i+1:nr_sites
            nr_jumps = sites.transitions(i,j) + sites.transitions(j,i);
            if nr_jumps > 0
                dist = sqrt(calc_dist_sqrd_frac(sites.frac_pos(:,i), sites.frac_pos(:,j), lattice)); 
                % For the jump diffusivity:
                jump_diff = jump_diff + nr_jumps*(dist^2);
                % For the histogram:
                dist_bin = ceil(dist/resolution);
                jumps_vs_dist(dist_bin) = jumps_vs_dist(dist_bin) + nr_jumps;
                if dist_bin > max_bin
                    max_bin = dist_bin;
                end
            end
        end
    end

    fprintf('Jump diffusion calculated assuming %f dimensional diffusion. \n', dimensions)
    jump_diff = (jump_diff*ang2m2)/(2*dimensions*nr_diffusing*sim_time); %In m^2/sec
    fprintf('Jump diffusivity (in m^2/sec): %e \n', jump_diff)
    
    %% Plot histogram of jumps vs distance of the jump
    if pics 
        figure
        distances = (0.5*resolution:resolution:(max_bin+0.5)*resolution);
        bar(distances, jumps_vs_dist(1:max_bin+1))
        title('Jumps vs. distance')
        xlabel('Distance (Angstrom)')
        ylabel('Nr. of jumps')
    end
end