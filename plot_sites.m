function plot_jump_paths(sites, lat, show_names)
% Plot the sites and the jumps between them
    nr_sites = size(sites.occupancy, 1); 
    figure;
    hold on
   
%% Plot the lattice:
    for i = 1:3 
        plot3([0, lat(i,1)], [0, lat(i,2)], [0, lat(i,3)], 'Color', 'black') %from the origin    
        for j = 1:3 % the combinations of two vectors
            if i ~= j
                plot3([lat(i,1) lat(i,1)+lat(j,1)], [lat(i,2) lat(i,2)+lat(j,2)], ...
                    [lat(i,3) lat(i,3)+lat(j,3)], 'Color', 'black')        
                for k = 1:3
                    if k ~= i && k ~= j % Combinations of three vectors:
                        plot3([lat(i,1)+lat(j,1), lat(i,1)+lat(j,1)+lat(k,1)], ...
                            [lat(i,2)+lat(j,2), lat(i,2)+lat(j,2)+lat(k,2)], ...
                            [lat(i,3)+lat(j,3), lat(i,3)+lat(j,3)+lat(k,3)], 'Color', 'black')
                    end
                end
            end            
        end
    end

%% Plot the sites (and site-names if show_names == True)
delta = 0.2; %Small displacement for better text visibility for the names
for i = 1:nr_sites
    point_size = 15*log(sites.occupancy(i)) + 1;
    if point_size > 1
        scatter3(sites.cart_pos(1,i), sites.cart_pos(2,i), sites.cart_pos(3,i), ...
        point_size, 'blue', 'filled');
        if show_names
            text(sites.cart_pos(1,i)+delta, sites.cart_pos(2,i)+delta, ...
                sites.cart_pos(3,i)+delta, sites.site_names(i) );          
        %else %show site_numbers
        %    text(sites.cart_pos(1,i)+delta, sites.cart_pos(2,i)+delta, ...
        %        sites.cart_pos(3,i)+delta, num2str(i) );  
        end
    end
end

% Plot the jumps as lines between the sites, and linewidth depending on the
% number of jumps
vec_cart = zeros(3,2);
vec_frac = zeros(3,2);
vec1 = zeros(3,2);
vec2 = zeros(3,2);

for i = 1:nr_sites-1
    for j = i+1:nr_sites
        linewidth = 0;
        if sites.transitions(i,j) + sites.transitions(j,i) > 0
            linewidth = 1 + log(sites.transitions(i,j) + sites.transitions(j,i));
        end
        if linewidth > 0 %plot the line
            color = 'red';
            %plot a line between pos(i) and pos(j), also take care of PBC
            for k = 1:3
                vec_cart(k,:) = [sites.cart_pos(k,i), sites.cart_pos(k,j)];
                vec_frac(k,:) = [sites.frac_pos(k,i), sites.frac_pos(k,j)];
            end
            pbc(1:3) = false ;
            nr_pbc = 0;
            for k = 1:3 %implement PBC
                if abs(sites.frac_pos(k,i) - sites.frac_pos(k,j)) > 0.5 %sites.positions(i,k) - sites.positions(j,k)) > 0.5
                    pbc(k) = true ;
                    if nr_pbc == 0
                        nr_pbc = 1;
                        [~, index_max] = max(vec_frac(k,:));
                        [~, index_min] = min(vec_frac(k,:));
                        vec1(k,:) = [max(vec_frac(k,:)), min(vec_frac(k,:))+1];
                        vec2(k,:) = [max(vec_frac(k,:))-1, min(vec_frac(k,:))];
                    else
                        [~, index_max2] = max(vec_frac(k,:));
                        if index_max2 == index_max %do the same as above
                            vec1(k,:) = [max(vec_frac(k,:)), min(vec_frac(k,:))+1];
                            vec2(k,:) = [max(vec_frac(k,:))-1, min(vec_frac(k,:))];    
                        else %do it the other way around
                            vec1(k,:) = [max(vec_frac(k,:))-1, min(vec_frac(k,:))];
                            vec2(k,:) = [max(vec_frac(k,:)), min(vec_frac(k,:))+1];                                 
                        end
                    end
                end
            end
            if pbc(1) || pbc(2) || pbc(3) %if PBC in one of the directions: plot vec1, vec2
                for k = 1:3
                   if not(pbc(k))                     
                       vec1(k,:) = [vec_frac(k,index_max), vec_frac(k,index_min)];
                       vec2(k,:) = [vec_frac(k,index_max), vec_frac(k,index_min)];
                   end 
                end            
                vec1a = frac_to_cart(vec1(:,1), lat); %Convert to cartesian coor
                vec1b = frac_to_cart(vec1(:,2), lat);
                vec2a = frac_to_cart(vec2(:,1), lat);
                vec2b = frac_to_cart(vec2(:,2), lat);
                plot3([vec1a(1) vec1b(1)], [vec1a(2) vec1b(2)], [vec1a(3) vec1b(3)], 'Color',color, 'LineWidth',linewidth);                  
                plot3([vec2a(1) vec2b(1)], [vec2a(2) vec2b(2)], [vec2a(3) vec2b(3)], 'Color',color, 'LineWidth',linewidth);                    
            else %plot vec
                plot3(vec_cart(1,:),vec_cart(2,:),vec_cart(3,:), 'Color',color, 'LineWidth',linewidth);
            end
        end
    end
end

% Make the plot look nice:
    title('Jumps between sites')
    xlabel('X (Angstrom)'); ylabel('Y (Angstrom)'); zlabel('Z (Angstrom)')
    axis equal
    min_coor = sum(lat.*(lat < 0)); % The minimum possible x,y,z coors
    max_coor = sum(lat.*(lat > 0)); % The maximum possible x,y,z coors
    axis([min_coor(1) max_coor(1) min_coor(2) max_coor(2) min_coor(3) max_coor(3)])
    view(-10, 10)
    hold off

end