function compare_sims()
% Compare the results (per property) from several different foldes and temperatures
% Combine upfolder, subs and temps to loop over all folders
% The results are saved in 'compare_file', from which the plots are made

    % Folders and subfolders to read:
    folder = 'Li3PS4/';
    subfolder = {'22Li/', '24Li/', '26Li/', 'Br22Li/', 'O24Li/'};
    temps = {'450K', '600K', '750K'};
    
    compare_file = [folder, 'sims_compare.mat'];
    % Load or construct the comparison file:
    if ~exist(compare_file, 'file')
        disp('sims_compare.mat not found, reading the data from other files')
        sims_comp = read_in_sims(folder, subfolder, temps);
        save(compare_file, 'sims_comp')
    else
        disp('sims_compare found')
        load(compare_file)
    end
    
    %Plot the results in the compare_file:
    plot_comparison(sims_comp, subfolder) 
end

function plot_comparison(sims_comp, subs)
% Plots all the properties from props_to_plot vs. temperature for the
% simulations given in sims_comp

    %Properties with 1 value per simulation:
    props_to_plot = {'vibration_amp', 'attempt_freq',  ...
        'tracer_diffusion', 'tracer_conductivity', ...  
        'total_occup', 'frac_collective', 'jump_diffusion', 'correlation_factor'};
    titles_of_plots = {'Vibration amplitude (Angstrom)', 'Attempt frequency (Hz)',  ...
       'Tracer diffusivity (m^2/sec)', 'Tracer conductivity (S/m)',...   
       'Known site occupation (%)', 'Collective jumps (%)', ...
       'Jump diffusivity (m^2/sec)', 'Correlation factor'};

   % Properties with a value per type of jump: 
    multi_props_to_plot  = {'e_act', 'rates'};
    multi_titles_of_plots = {'Activation energy (eV)', 'Jump rate (Hz)'};
   
    linestyles = {'-o', '-^', '-*', '-p', '-+'};
    pointstyles = {'+', 'o', '*', '+', 'o', '*', '+', 'o', '*', '+', 'o', '*', '+', 'o', '*'}
    
    sims = fieldnames(sims_comp);
    for i = 1:numel(sims) 
        temp = strsplit(sims{i}, '_');
        all_names{i} = temp{2}; % !! Assuming the 'material names' are in temp{2} !!
    end
    % Find the same names to plot those with a line between them:
    names = unique(all_names);
    %% Plot properties with 1 value per simulation versus Temperature
    for a = 1:numel(props_to_plot)
        figure()
        hold on
        % Labels:
        xlabel('Temperature (K)')
        ylabel(titles_of_plots{a})        
        for i = 1:numel(names)
            counter = 1;
            temp_x = 0;
            temp_y = 0;
            for j = 1:numel(all_names)
                if strcmp(names{i}, all_names{j}) %same names
                    temp_x(counter) = sims_comp.(sims{j}).temperature;
                    temp_y(counter) = sims_comp.(sims{j}).(props_to_plot{a});
                    counter = counter + 1;
                end
            end
            plot(temp_x, temp_y, linestyles{i}, 'LineWidth', 2.0, 'MarkerSize', 10.0)
        end       
        legend(names)
        % For checking if the legend is correct:        
        %names 
        % Log-scale is better in some cases:
        if strcmp(props_to_plot{a},'tracer_diffusion') || strcmp(props_to_plot{a}, 'tracer_conductivity') ...
            || strcmp(props_to_plot{a}, 'jump_diffusion')    
            set(gca, 'YScale', 'log')
        end
        hold off
    end
    
    %% Properties with a value per type of jump, plot versus jump name 
    for a = 1:numel(multi_props_to_plot)
        figure()
        ax = gca;
        hold on
        ylabel(multi_titles_of_plots{a})
        temp_x = [1:1:numel(sims_comp.(sims{1}).jump_names)];
        for i = 1:numel(sims)
            temp_y = sims_comp.(sims{i}).(multi_props_to_plot{a})(:,1)';
            plot(temp_x, temp_y, pointstyles{i}, 'LineWidth', 2.0, 'MarkerSize', 10.0)
        end       
        legend(subs)
        ax.XTick = temp_x; 
        % !! Assuming the same jump names in all simulations being compared !!
        ax.XTickLabels = sims_comp.(sims{1}).jump_names
        ax.XTickLabelRotation = 90;   
        grid('on')
        % Log scale is better in some cases:
        if strcmp(multi_props_to_plot{a}, 'rates')
            set(gca, 'YScale', 'log')
        end
        hold off
    end
end
    
function sims_comp = read_in_sims(upfolder, subs, temps)
% Combine upfolder, subs and temps to loop over all folders
% and read in certain data of all the given simulations to compare them
% easily using plot_comparisons
    for i = 1:size(subs,2)
        for j = 1:size(temps,2)
            folder = [upfolder, subs{i}, temps{j}];
            sim_data_file = [folder, '/simulation_data.mat'];
            sites_file = [folder, '/sites.mat'];            
            if exist(sim_data_file, 'file') &&  exist(sites_file, 'file')
                fprintf('Loading simulation_data.mat and sites.mat from %s \n', folder)
                load(sim_data_file) 
                load(sites_file)
                % The info from this simulation:
                info = read_sim_info(sim_data, sites);
                % The 'name' of the simulation:
                sim_name = strrep(folder, '/', '_');
                sims_comp.(sim_name) = info;
            elseif exist(sim_data_file, 'file')
                fprintf('sites.mat NOT found in given folder: %s \n', folder)                   
            else
                fprintf('sim_data.mat and sites_sites.mat NOT found in given folder: %s \n', folder)
            end
        end %j
    end
    
end

function info = read_sim_info(sim_data, sites)
% The information to be compared for the simulations:
% All the info wanted from sim_data:
    info.temperature = sim_data.temperature;
    info.attempt_freq = sim_data.attempt_freq;
    info.vibration_amp = sim_data.vibration_amp;
    info.tracer_diffusion = sim_data.tracer_diffusion;
    info.tracer_conductivity = sim_data.tracer_conductivity;
    
% All the info wanted from sites:    
    %Fraction of time the atom is at known locations:
    info.total_occup = 100*sum(sites.occupancy)/numel(sites.atoms); 
    info.site_occup = 100.*sites.sites_occup;
    info.atom_locations = 100.*sites.atom_locations;
    info.site_names = sites.stable_names;
    % Jump names: 
    info.jump_names = sites.jump_names;
    % Activation energy
    info.e_act = sites.e_act;
    % Jump rates:
    info.rates = sites.rates;   
    info.jump_diffusion = sites.jump_diffusivity;
    info.correlation_factor = sites.correlation_factor;
    % Fraction of collective jumps
    info.frac_collective = 100*(1.0-sites.solo_frac);
end