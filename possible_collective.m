function [collective, coll_jumps, coll_matrix, multi_coll, uncoll_count] = ... 
    possible_collective(sites, sim_data, coll_dist)
% Determine the amount of jumps (and which) could show collective behaviour
    disp('Determining collective jumps')
    % The number of steps which would still 'mean' correlation:
    coll_steps = ceil((1.0/sim_data.attempt_freq)/sim_data.time_step);

    coll_dist_sqrd = coll_dist*coll_dist; %for easy comparison
    lattice = sim_data.lattice;

    % The times at which jumps occur:
    if sum(sum(sites.transitions)) > 0
        [times, orig] = sort(sites.all_trans(:,5));   
        nr_jumps = length(times);
    else
        nr_jumps = 0;
    end
    
    % define some stuff
    coll_count = 0;
    uncoll_count = 0;
    collective(1, :) = [0, 0];     
    coll_jumps = {'', ''};
    if nr_jumps > 1 % For only 1 jump no correlations can occur
        for i = 1:nr_jumps-1 % -1 because no need to check the last one
            j = i + 1;
            while j <= nr_jumps && (times(j)-times(i) <= coll_steps) %time correlation is ok..
                if sites.all_trans(orig(i),1) ~= sites.all_trans(orig(j),1) ... %different atoms
                    % The site numbers of start and end sites 
                    start_site_i = sites.all_trans(orig(i),2); %start sites
                    start_site_j = sites.all_trans(orig(j),2);
                    end_site_i = sites.all_trans(orig(i),3); %end sites
                    end_site_j = sites.all_trans(orig(j),3);                    
                    % Check distances between sites (all combinations: ss, se (= start-end),
                    % es, ee)
                    ss_dist_sqrd = calc_dist_sqrd_frac(sites.frac_pos(:,start_site_i), sites.frac_pos(:,start_site_j), lattice);
                    se_dist_sqrd = calc_dist_sqrd_frac(sites.frac_pos(:,start_site_i), sites.frac_pos(:,end_site_j), lattice);
                    es_dist_sqrd = calc_dist_sqrd_frac(sites.frac_pos(:,end_site_i), sites.frac_pos(:,start_site_j), lattice);
                    ee_dist_sqrd = calc_dist_sqrd_frac(sites.frac_pos(:,end_site_i), sites.frac_pos(:,end_site_j), lattice);
                    if ss_dist_sqrd <= coll_dist_sqrd || se_dist_sqrd <= coll_dist_sqrd || ... 
                        es_dist_sqrd <= coll_dist_sqrd || ee_dist_sqrd <= coll_dist_sqrd
                    %fprintf('Possibly correlated jump detected at timestep %d, atoms involved are: %d and %d \n',...
                    %    times(i), sites.all_trans(orig(i),1), sites.all_trans(orig(j),1))
                    % remember the numbers in all_trans of these two jumps
                        coll_count = coll_count + 1;
                        collective(coll_count, :) = [orig(i), orig(j)] ;
                        coll_jumps(coll_count, :) = {[sites.site_names{start_site_i}, '_to_', sites.site_names{end_site_i}], ... 
                            [sites.site_names{start_site_j}, '_to_', sites.site_names{end_site_j}]}; 
                    else
                        uncoll_count = uncoll_count + 1;
                    end
                else
                    uncoll_count = uncoll_count + 1;
                end
                j = j + 1;
            end
        end
    end
    fprintf('Total number of jumps: %d \n', nr_jumps)
    fprintf('Number of possibly collective jumps found: %d \n', coll_count)
    fprintf('Number of solo jumps found: %d \n', uncoll_count)

    % Analyse what types of jumps are collective and make a matrix of all jump combinations.
    coll_matrix = zeros(size(sites.jump_names,1));
    if coll_count > 0
        for a = 1:coll_count
            i = 0;
            j = 0;
            for b = 1:size(sites.jump_names,1)
                if strcmp(sites.jump_names{b}, coll_jumps{a,1})
                    i = b; %the first jump name
                end
                if strcmp(sites.jump_names{b}, coll_jumps{a,2})
                    j = b; %the second jump name
                end
            end
            coll_matrix(i,j) = coll_matrix(i,j) + 1;
        end
    end
  
    multi_coll = 0;
    if coll_count > 1 % Check if there are multi-collective jumps...
        uni = unique(collective);
        if numel(uni) ~= numel(collective) % Some jumps are 'multi' correlated
            fprintf('Nr. of multiply collective jumps found: %d \n', numel(collective) - numel(uni))
            % Find which ones occur multiple times, and with which others....
            multi_coll = zeros(numel(collective) - numel(uni),1);
            count = 1;
            sorted = sort([collective(:,1); collective(:,2)]);
            for i = 2:numel(sorted)
                if sorted(i) == sorted(i-1)
                    multi_coll(count) = sorted(i);
                    count = count + 1;
             %       disp('Non-unique jump found!')
                end
            end
        end
    end
    
    
end