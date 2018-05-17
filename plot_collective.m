function plot_collective(sites)
%% Plot a histogram of the number of jumps vs. timestep
    centers = [250:500:size(sites.atoms(:,1))];
    times = sort(sites.all_trans(:,5));   
    figure
    histogram(times, centers) 
    title('Histogram of jumps vs. time')
    xlabel('Time (steps)')
    ylabel('Nr. of jumps')    
    
%% Plot which types of jumps are 'collective'    
    figure
    image(sites.coll_matrix, 'CDataMapping','scaled') 
    ax = gca;
    colorbar
    title('Number of cooperative jumps per jump-type combination')
    % Put all the jump names on the x- and y-axis
    ax.XTick = 1:1:size(sites.jump_names,1); 
    ax.YTick = 1:1:size(sites.jump_names,1);
    ax.XTickLabels = sites.jump_names;
    ax.YTickLabels = sites.jump_names;
    ax.XTickLabelRotation = 90;  
end