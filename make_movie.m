function make_movie(sites, sim_data, start_end, movie_file, step_frame)
    % Make a movie showing how the atoms jump from site to site
    % LATTICE NOT TAKEN INTO ACCOUNT CORRECTLY YET!
    %tic
    
    nr_atoms = sim_data.nr_diffusing;
    F = VideoWriter(movie_file, 'MPEG-4');
    open(F);
    
    % FOR SOME STRANGE REASON SHOWING THE MOVIE TO THE SCREEN IS ~4 TIMES
    % FASTER AS NOT SHOWING IT...
    %figurehandle = figure('Visible','off'); %To not show the plots while making a movie
    figurehandle = figure(); %Show the plot
    % The look of the figure:    
    hold on;
    xlabel('a'); 
    ylabel('b'); 
    zlabel('c');
    axis equal;
    axis([0 1 0 1 0 1]);
    ax = gca ;
    ax.XTick = [0 0.25 0.5 0.75 1]; 
    ax.YTick = [0 0.25 0.5 0.75 1];
    ax.ZTick = [0 0.25 0.5 0.75 1];
    view(-10,10); % vary with loop if you want plot to rotate during video
    grid on;
    set(gcf); 
    
% Initialise some stuff:
    color_vec = hsv(nr_atoms); % Generate vector with colors for each atom
    line_color = 'red'; %Color of the connecting lines
    linewidth = 3.0; % Size of the connecting lines
    pos_size = 200; % Size of the positions
    atom_size = 170; % Size of the atoms
    
    disp('Making a movie of the jumps, will take some time (depending on the number of jumps)...')
    fprintf('Number of time steps per frame: %d \n', step_frame)
    
    % Get frames at which a jump occurs and where it is in all_trans   
    [all_jump_times, trans_nr] = sort(sites.all_trans(:,5));
    all_jumps = size(trans_nr, 1);
    start_jump = 0;
    end_jump = 0;
    for i = 1:all_jumps
        if all_jump_times(i) >  start_end(1) && start_jump == 0  
            start_jump = i;
        elseif all_jump_times(i) > start_end(2) && end_jump == 0
            end_jump = i;
        end
    end
    
    if start_jump == 0 
        disp('ERROR! No jumps between the given time steps, not making a movie!')
        return
    elseif end_jump == 0
        end_jump = all_jumps; 
    end
    jump_times = all_jump_times(start_jump:end_jump);
    nr_jumps = size(jump_times, 1);
    
    nr_stable = size(sites.frac_pos,2);
    lines1 = plot3([0 0], [0 0], [0 0], 'Color', line_color, 'LineWidth', 1E-10);
    lines2 = plot3([0 0], [0 0], [0 0], 'Color', line_color, 'LineWidth', 1E-10);
    
    % Draw first frame:
    % Draw the positions
    scatter3(sites.frac_pos(1,:),sites.frac_pos(2,:),sites.frac_pos(3,:), pos_size,'MarkerEdgeColor','black','LineWidth',1);
    % The atoms:
    atoms = scatter3(zeros(nr_atoms,1),zeros(nr_atoms,1),zeros(nr_atoms,1), atom_size, 'CData', color_vec, 'MarkerFaceColor', 'flat');    
    for j = 1:nr_atoms   
        site_nr = sites.atoms(start_end(1),j); %The position at the chosen start time
        if site_nr > 0 && site_nr <= nr_stable %Draw the starting position
            atoms.XData(j) = sites.frac_pos(1,site_nr);
            atoms.YData(j) = sites.frac_pos(2,site_nr);
            atoms.ZData(j) = sites.frac_pos(3,site_nr);
        else %Find the first position this atom occupies:
            t = start_end(1);
            while(site_nr == 0 || site_nr > nr_stable)
                t = t + 1;
                site_nr = sites.atoms(t,j);
            end
            atoms.XData(j) = sites.frac_pos(1,site_nr);
            atoms.YData(j) = sites.frac_pos(2,site_nr);
            atoms.ZData(j) = sites.frac_pos(3,site_nr);
        end
    end
    writeVideo(F,getframe(figurehandle));
    
    % Make the frames
    jump_count = 1;
    disp_count = 0;
    for i = start_end(1):start_end(2)
        while jump_count <= nr_jumps && jump_times(jump_count) == i % Allows for multiple jumps at the same time-step
                % the jump:
                jump_atom = sites.all_trans(trans_nr(jump_count),1);
                prev_site = sites.all_trans(trans_nr(jump_count),2);
                site_nr = sites.all_trans(trans_nr(jump_count),3);
                prev_pos = [sites.frac_pos(1,prev_site) sites.frac_pos(2,prev_site) sites.frac_pos(3,prev_site)];             
                % Draw the atom at the new positions
                atoms.XData(jump_atom) = sites.frac_pos(1,site_nr);
                atoms.YData(jump_atom) = sites.frac_pos(2,site_nr);
                atoms.ZData(jump_atom) = sites.frac_pos(3,site_nr);
                % Draw a line between the begin and end positions
                new_pos = [atoms.XData(jump_atom) atoms.YData(jump_atom) atoms.ZData(jump_atom)]; 
                [vec1, vec2] = line_pbc(prev_pos, new_pos);
                if any(vec2) %vec2 contains a vector, so plot both
                    lines1(jump_count) = plot3(vec1(1,:),vec1(2,:),vec1(3,:), 'Color', line_color, 'LineWidth', linewidth);
                    lines2(jump_count) = plot3(vec2(1,:),vec2(2,:),vec2(3,:), 'Color', line_color, 'LineWidth', linewidth);                 
                    %lines(length(lines)+1) = plot3(vec1(1,:),vec1(2,:),vec1(3,:), 'Color', line_color, 'LineWidth', linewidth);
                    %lines(length(lines)+1) = plot3(vec2(1,:),vec2(2,:),vec2(3,:), 'Color', line_color, 'LineWidth', linewidth);
                else % vec2 is empty
                    %lines(length(lines)+1) = plot3(vec1(1,:),vec1(2,:),vec1(3,:), 'Color', line_color, 'LineWidth', linewidth);
                    lines1(jump_count) = plot3(vec1(1,:),vec1(2,:),vec1(3,:), 'Color', line_color, 'LineWidth', linewidth);
                    % QUICK FIX:
                    lines2(jump_count) = plot3([0 0],[0 0],[0 0], 'Color', line_color, 'LineWidth', 1E-10);  
                end
            jump_count = jump_count + 1;     
        end
        if mod(i, step_frame) == 0 
            % Every frame reduce the linewidth of all the lines by 'disappear_speed' to have them
            % slowly disappear. 0.0625: After 2 sec (with a starting linewidth of 3 and 24 frames/sec)
            disappear_speed = 0.0625 - 1E-10/48; %-1E-10 to avoid 'line width smaller as zero' errors
            for j = disp_count+1:jump_count-1
                if lines1(j).LineWidth > 1E-9
                    lines1(j).LineWidth = lines1(j).LineWidth - disappear_speed;
                else % This line has 'disappeared' --> no need to check anymore, update disp_count
                    disp_count = disp_count + 1;
                end
                if lines2(j).LineWidth > 1E-9
                    lines2(j).LineWidth = lines2(j).LineWidth - disappear_speed;
                end
            end           
            title(['Timestep ', int2str(i)])
            writeVideo(F,getframe(figurehandle));
        end
    end
    %Done...
    close(F)
    %toc
    disp('Movie done')
    
end

function [vec1, vec2] = line_pbc(pos1, pos2)
%plot a line between pos1 and pos2, also take care of PBC            
    vec(1,:) = [pos1(1),pos2(1)];
    vec(2,:) = [pos1(2),pos2(2)];
    vec(3,:) = [pos1(3),pos2(3)];
    pbc = zeros(3,1); %= False
    vec1 = zeros(3,2);
    vec2 = zeros(3,2);
    
    nr_pbc = 0;
    for k = 1:3 %implement PBC
        if abs(vec(k,1) - vec(k,2)) > 0.5
            pbc(k) = true ;
            if nr_pbc == 0
                nr_pbc = 1;
                [~, index_max] = max(vec(k,:));
                [~, index_min] = min(vec(k,:));
                vec1(k,:) = [max(vec(k,:)), min(vec(k,:))+1];
                vec2(k,:) = [max(vec(k,:))-1, min(vec(k,:))];
            else
                [~, index_max2] = max(vec(k,:));
                if index_max2 == index_max %do the same as above
                    vec1(k,:) = [max(vec(k,:)), min(vec(k,:))+1];
                    vec2(k,:) = [max(vec(k,:))-1, min(vec(k,:))];    
                else %do it the other way around
                    vec1(k,:) = [max(vec(k,:))-1, min(vec(k,:))];
                    vec2(k,:) = [max(vec(k,:)), min(vec(k,:))+1];                                 
                end
            end
        end
    end
% if PBC in one of the directions make vec1 and vec2    
    if pbc(1) || pbc(2) || pbc(3) 
        for k = 1:3
            if not(pbc(k))                     
                vec1(k,:) = [vec(k,index_max), vec(k,index_min)];
                vec2(k,:) = [vec(k,index_max), vec(k,index_min)];
            end 
        end                  
    else %No pbc --> vec1 = vec
        vec1 = vec;
    end
    
end

