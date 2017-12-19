function plot_rdfs(rdf)
% Plot all the rdf's
    nr_rdfs = numel(rdf.rdf_names);
    
    for i = 1:nr_rdfs
        if rdf.integrated(1,end,i) > 0
            % Density vs. distance
            figure
            plot(rdf.distributions(:,:,i)')
            title(rdf.rdf_names{i}) 
            xticklabels(0:rdf.max_dist)
            xlabel('Distance (Angstrom)')
            ylabel('Density (a.u)')
            legend(rdf.elements) 
            grid on
            
            % The INTEGRATED density vs. distance:
            %figure
            %plot(rdf.integrated(:,:,i)')
            %title(rdf.rdf_names{i})   
            %xticklabels(0:rdf.max_dist)
            %xlabel('Distance (Angstrom)')
            %ylabel('Density (a.u)')
            %legend(rdf.elements)
            %grid on
        end
    end