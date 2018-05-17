function plots_from_displacement(sim_data, sites, resolution)
% Plots based on information from the displacement (mainly)
%% Define some stuff:
%% At which nr. the diffusing atoms start and end:
    start_diff = sim_data.start_diff_elem;
    end_diff = sim_data.end_diff_elem;

    %% Plot displacement for each diffusing atom:
    figure
    hold on
    for i = start_diff:end_diff 
        plot(sim_data.displacement(i,:))
    end
    hold off
    title('Displacement of diffusing element')
    xlabel('Time step')
    ylabel('Displacement (Angstrom)')
    
    %% Plot histogram of final displacement of diffusing atoms:
    figure
    hist(sim_data.displacement(start_diff:end_diff,sim_data.nr_steps)')
    title('Histogram of displacement of diffusing element')
    xlabel('Displacement (Angstrom)')
    ylabel('Nr. of atoms')
    
    %% Plot displacement per element: 
    counter = 1;
    disp_elem = zeros(sim_data.nr_elements, sim_data.nr_steps);
    figure
    hold on
    for i = 1:sim_data.nr_elements
        for j = counter:(counter+sim_data.nr_per_element(i)-1)
            disp_elem(i,:) = disp_elem(i,:) + sim_data.displacement(j,:);
        end
        counter = counter + sim_data.nr_per_element(i);
        disp_elem(i,:) = disp_elem(i,:)./sim_data.nr_per_element(i);
        plot(disp_elem(i,:))
    end
    hold off
    title('Displacement per element')
    xlabel('Time step')
    ylabel('Displacement (Angstrom)')
    legend(sim_data.elements) 
    
    %% Density plot of diffusing element
    lat = sim_data.lattice;
    length = zeros(3,1);
    % Find minimal and maximal x,y,z values in cart_pos
    min_coor = sum(lat.*(lat < 0)); % The minimum possible x,y,z coors
    max_coor = sum(lat.*(lat > 0)); % The maximum possible x,y,z coors
    for i = 1:3
        length(i) = ceil((max_coor(i) - min_coor(i))/resolution);
    end
    
    % Plot the lattice:
    figure
    hold on
    origin = [-min_coor(1); -min_coor(2); -min_coor(3)]; %shift origin to easily calculate the density postions
    lat = lat./resolution;
    origin = origin./resolution;
    lat_orig = zeros(3);
    for i = 1:3
        for j = 1:3
            lat_orig(i,j) = lat(i,j) + origin(j);
        end
    end

    for i = 1:3 %from the origin
        plot3([origin(1), lat_orig(i,1)], [origin(2), lat_orig(i,2)], ...
            [origin(3), lat_orig(i,3)], 'Color', 'black')     
        for j = 1:3 % the combinations of two vectors
            if i ~= j
                plot3([lat_orig(i,1) lat_orig(i,1)+lat(j,1)], [lat_orig(i,2) lat_orig(i,2)+lat(j,2)], ...
                    [lat_orig(i,3) lat_orig(i,3)+lat(j,3)], 'Color', 'black')        
                for k = 1:3
                    if k ~= i && k ~= j 
                        plot3([lat_orig(i,1)+lat(j,1), origin(1)+lat(i,1)+lat(j,1)+lat(k,1)], ...
                            [lat_orig(i,2)+lat(j,2), origin(2)+lat(i,2)+lat(j,2)+lat(k,2)], ...
                            [lat_orig(i,3)+lat(j,3), origin(3)+lat(i,3)+lat(j,3)+lat(k,3)], 'Color', 'black')
                    end
                end
            end            
        end
    end
    
    box = zeros(length(1), length(2), length(3));
    % Now fill the box with densities:
    pos = zeros(3,1);
    for atom = start_diff:end_diff 
        for time = 1:sim_data.nr_steps
            for i = 1:3
            % Calculate the position in the box
                pos(i) = 1 + floor((sim_data.cart_pos(i,atom,time) - min_coor(i))/resolution);
            end
            box(pos(1),pos(2),pos(3)) = box(pos(1),pos(2),pos(3)) + 1;
        end
    end
    % Show it more nicely by smoothing the density:
    box_smooth = smooth3(box, 'gaussian');
    
    maxval = max(max(max(box_smooth))); % the maximum value of the density
    % Define the color, isovalue & alpha for the density plot: 
    colors = ['r' 'y' 'c']; % red, yellow, cyan
    iso = [0.25 0.10 0.02]; % isosurface values (with respect to the maximum value)
    alpha = [0.5 0.3 0.1]; % the transparancy;
    % For some reason the x and y axis are permuted by isosurface, fix this by permuting
    % them again:
    box_smooth = permute(box_smooth, [2, 1, 3]);
    for i = 1:size(iso,2)
        isoval = iso(i)*maxval;
        surf=isosurface(box_smooth,isoval); 
        p1 = patch(surf);
        isonormals(box_smooth,p1);
        % set the color, mesh and transparency level of the surface
        set(p1,'FaceColor',colors(i),'EdgeColor','none','FaceAlpha',alpha(i)); 
    end  
    
    %% Plot the site positions on top of the density
    delta = 0.2; %Small displacement for better text visibility
    nr_sites = size(sites.transitions,1);
    for i = 1:nr_sites
        point_size = 100;
        % Plot site positions:
        scatter3(origin(1) + sites.cart_pos(1,i)/resolution, ...
            origin(2) + sites.cart_pos(2,i)/resolution, ...
            origin(3) + sites.cart_pos(3,i)/resolution, point_size, 'blue', 'filled');
        % Add name to the positions:
        text(origin(1)+(sites.cart_pos(1,i)+delta)/resolution, ...
             origin(2)+(sites.cart_pos(2,i)+delta)/resolution, ...
             origin(3)+(sites.cart_pos(3,i)+delta)/resolution, sites.site_names(i) );  
        % OR show site_numbers instead of names
        %text(origin(1)+(sites.cart_pos(1,i)+delta)/resolution, ...
        %    origin(2)+(sites.cart_pos(2,i)+delta)/resolution, ...
        %    origin(3)+(sites.cart_pos(3,i)+delta)/resolution, num2str(i) );  
    end    
    
 %% Add the STARTING positions of the non-diffusing atoms, if necessary
    % DIFFERENT COLORS FOR DIFFERENT ELEMENTS?
%     time = 1; % The starting position, just to give an estimate
%     for i = 1:sim_data.nr_atoms
%         if i < sim_data.start_diff_elem || i > sim_data.end_diff_elem
%              % Plot atom positions:
%             point_size = 250;
%             if strcmp('S', sim_data.atom_element(i)) %Sulfur
%                 color = 'green';
%             elseif strcmp('O', sim_data.atom_element(i)) %Oxygen
%                 color = 'magenta';
%             elseif strcmp('P', sim_data.atom_element(i)) %Phosphor
%                 color = 'black';
%             else
%                 color = 'red';
%             end    
%             % The position
%             scatter3(origin(1) + sim_data.cart_pos(1,i,time)/resolution, ...
%                 origin(2) + sim_data.cart_pos(2,i,time)/resolution, ...
%                 origin(3) + sim_data.cart_pos(3,i,time)/resolution, point_size, color, 'filled');
%             % The element, a bit much to always show
%             text(origin(1) + (sim_data.cart_pos(1,i,time)+delta)/resolution, ...
%                 origin(2) + (sim_data.cart_pos(2,i,time)+delta)/resolution, ...
%                 origin(3) + (sim_data.cart_pos(3,i,time)+delta)/resolution, sim_data.atom_element(i));   
%         end
%     end
 %% Add title, labels and ticks:
    hold off
    title('Density of diffusing element')
    axis equal
    axis([0 length(1) 0 length(2) 0 length(3)])
    xlabel('X (Angstrom)'); ylabel('Y (Angstrom)'); zlabel('Z (Angstrom)')
    ax = gca ;
    % Add tick labels every 'angstrom_spacing' Angstrom:
    angstrom_spacing = 2;
    label_spacing = round(angstrom_spacing/resolution); %Nr. of steps between ticks
    ax.XTick = 0:label_spacing:length(1);
    ax.YTick = 0:label_spacing:length(2);
    ax.ZTick = 0:label_spacing:length(3);
    ax.XTickLabels = round(min_coor(1)):angstrom_spacing:round(max_coor(1));
    ax.YTickLabels = round(min_coor(2)):angstrom_spacing:round(max_coor(2));
    ax.ZTickLabels = round(min_coor(3)):angstrom_spacing:round(max_coor(3));
    view(-10, 10)
    grid on

end    