function [names, pos, supercell] = known_materials(material)
% The positions of diffusing atoms, and their names, in different materials
    if strcmp(material, 'argyrodite')
       supercell = [1 1 1]; 
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3))    
       [names, pos] = argyrodite(supercell);
       
    elseif strcmp(material, 'latp')
       supercell = [1 1 1]; 
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3)) 
       [names, pos] = latp(supercell);
       
    elseif strcmp(material, 'na3ps4')
       supercell = [2 2 2];
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3)) 
       [names, pos] = na3ps4(supercell);
       
    elseif strcmp(material, 'lisnps')
       supercell = [1 1 1]; 
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3))         
       [names, pos] = lisnps(supercell);
       
    elseif strcmp(material, 'li3ps4_beta')
       supercell = [1 1 2]; 
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3))                 
       [names, pos] = li3ps4_beta(supercell);
       
    elseif strcmp(material, 'mno2_lambda')
       supercell = [1 1 1]; 
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3))    
       [names, pos] = mno2_lambda(supercell);             
       
    elseif strcmp(material, 'lagp')
       supercell = [1 2 2]; 
       fprintf('Material is: %s \n', material)         
       fprintf('Assuming a %d x %d x %d supercell \n', supercell(1), supercell(2), supercell(3))    
       [names, pos] = lagp(supercell);  
       
    else
       disp('ERROR! Material not recognised, add it to known_materials.m') 
       disp('ERROR! Also check if you have given the material name in all lower case correctly!')
    end
end

%% LAGP
function [names, pos] = lagp(supercell)
% Space group no. 167 (R-3c) h-axes
%The symmetries applicable to each site:
	sym = [0 0 0; 2/3 1/3 1/3; 1/3 2/3 2/3];
% 36f-sites
    x = 0.025;
    y = 0.025;
    z = 0.0;
    pos_sym = [ x,y,z;	-y,x-y,z;	-x+y,-x,z;	y,x,-z+1/2; ...
        x-y,-y,-z+1/2;	-x,-x+y,-z+1/2;	-x,-y,-z;	y,-x+y,-z; ...
    x-y,x,-z;	-y,-x,z+1/2;	-x+y,y,z+1/2;	x,x-y,z+1/2 ];
    names_sym = {'36f', '36f', '36f','36f', '36f', '36f','36f', '36f', '36f','36f', '36f', '36f'};

    %Construct all sites:
    [names, pos] = construct(sym', pos_sym', names_sym, supercell);
end

%% Na3PS4
function [names, pos] = na3ps4(supercell)
% Space group no. 217 (I-43m)
%The symmetries applicable to each site:
	sym = [0 0 0; 0.5 0.5 0.5];
% 6b-sites
    x = 0.5;
    y = 0.0;
    pos_sym = [ y x x; x y x; x x y ];
    names_sym = {'6b', '6b', '6b'};

    %Construct all sites:
    [names, pos] = construct(sym', pos_sym', names_sym, supercell);
end

%% NaMnO2_lambda
function [names,pos] = mno2_lambda(supercell)
% Spacegroup 227, Fd-3m, origin choice 2s
    disp('ORIGIN CHOICE 2 for NaMnO2 lambda')
%The symmetries applicable to each site:
    sym = [0 0 0; 0 0.5 0.5; 0.5 0 0.5; 0.5 0.5 0];
    % 16c Na-site
    pos_16c = [0 0 0; 0.75 0.25 0.5; 0.25 0.5 0.75; 0.5 0.75 0.25];
    names_16c = {'16c', '16c', '16c', '16c'}; 
    % 8a Na-site
    pos_8a = [0.125 0.125 0.125; 0.875 0.375 0.375];
    names_8a = {'8a', '8a'};    
   
  %All the positions in this material             
    pos_sym = [pos_16c(:,:); pos_8a(:,:)];
    names_sym = {names_16c{:}, names_8a{:} };   
  % Construct all:  
   [names, pos] = construct(sym', pos_sym', names_sym, supercell);
end

%% Li10SnP2S12
function [names,pos] = lisnps(supercell)
% Spacegroup 137, P4_2/nmc
    disp('ORIGIN CHOICE 2 for Li10SnP2S12!')
%The symmetries applicable to each site:
    sym = [0 0 0];
    % Li1-sites, 16h
    x = 0.514; y = 0.01; z = 0.072;
    pos_li1 = [x, y, z; (-x)+0.5, (-y)+0.5, z; (-y)+0.5, x, z+0.5; y, (-x)+0.5, z+0.5; ... 
        -x, y+0.5, -z; x+0.5, -y, -z; y+0.5, x+0.5, (-z)+0.5; -y, -x, (-z)+0.5; ...
        -x, -y, -z; x+0.5, y+0.5, -z; y+0.5, -x, (-z)+0.5; -y, x+0.5, (-z)+0.5; ...
        x, (-y)+0.5, z; (-x)+0.5, y, z; (-y)+0.5, (-x)+0.5, z+0.5; y, x, z+0.5];
    names_li1 = {'Li1', 'Li1', 'Li1', 'Li1', 'Li1', 'Li1', 'Li1', 'Li1', ...
        'Li1', 'Li1', 'Li1', 'Li1', 'Li1', 'Li1', 'Li1', 'Li1'}; 
    % Li2-sites, 16h
    x = 0.027; y = 0.007; z = 0.215;
    pos_li2 = [x, y, z; (-x)+0.5, (-y)+0.5, z; (-y)+0.5, x, z+0.5; y, (-x)+0.5, z+0.5; ... 
        -x, y+0.5, -z; x+0.5, -y, -z; y+0.5, x+0.5, (-z)+0.5; -y, -x, (-z)+0.5; ...
        -x, -y, -z; x+0.5, y+0.5, -z; y+0.5, -x, (-z)+0.5; -y, x+0.5, (-z)+0.5; ...
        x, (-y)+0.5, z; (-x)+0.5, y, z; (-y)+0.5, (-x)+0.5, z+0.5; y, x, z+0.5];
    names_li2 = {'Li2', 'Li2', 'Li2', 'Li2', 'Li2', 'Li2', 'Li2', 'Li2', ...
        'Li2', 'Li2', 'Li2', 'Li2', 'Li2', 'Li2', 'Li2', 'Li2'};    
    % Li3-sites, 4d
    z = 0.1937;
    pos_li3 = [0.25, 0.25, z; 0.25, 0.25, z+0.5; 0.75, 0.75, -z; 0.75, 0.75, (-z)+0.5];
    names_li3 = {'Li3', 'Li3', 'Li3', 'Li3'};     
    % Li4-sites, 4c
    z = 0.0042;
    pos_li4 = [0.75, 0.25, z; 0.25, 0.75, z+0.5; 0.25, 0.75, -z; 0.75, 0.25, (-z)+0.5];
    names_li4 = {'Li4', 'Li4', 'Li4', 'Li4'}; 
  %All the positions in this material             
    pos_sym = [pos_li1(:,:); pos_li2(:,:); pos_li3(:,:); pos_li4(:,:)];
    names_sym = {names_li1{:}, names_li2{:}, names_li3{:}, names_li4{:} };   
  % Construct all:  
   [names, pos] = construct(sym', pos_sym', names_sym, supercell);
end

%% Argyrodites:
function [names, pos] = argyrodite(supercell)
%The symmetries applicable to each site:
    sym = [0 0 0; 0 1/2 1/2; 1/2 0 1/2; 1/2 1/2 0];
%Only 48h positions:
    x = 0.183; 
    %y = x; 
    z = 0.024;
    pos_sym = [x x z; -x -x z; -x x -z; x -x -z; z x x; z -x -x; -z -x x; -z x -x; ... 
        x z x; -x z -x; x -z -x; -x -z x];
    names_sym = {'48h', '48h', '48h', '48h', '48h', '48h', '48h', '48h', '48h', ...
         '48h', '48h', '48h'};
    %Construct all sites:
    [names, pos] = construct(sym', pos_sym', names_sym, supercell);
end

%% LATP
function [names, pos] = latp(supercell)
%The symmetries applicable to each site, spacegroup 167 (R-3c)
    sym = [0 0 0; 2/3 1/3 1/3; 1/3 2/3 2/3];
%6b positions
    pos_6b = [0 0 0; 0 0 0.5];
    names_6b = {'6b', '6b'};
%18e positions (in between 36f's, roughly estimated, might need to be optimised)
    x = 0.66;
    pos_18e = [x, 0, 0.25; 0, x, 0.25; -x, -x, 0.25; -x, 0, 0.75; 0, -x, 0.75; x, x, 0.75];
    names_18e = {'18e', '18e', '18e', '18e', '18e', '18e'};
%36f positions, but it seems unrealistic to have them so close (within 0.5 Angstrom) together
%    x = 0.07; y = 0.34; z = 0.07;
%    pos_36f = [x, y, z; -y, x-y, z; (-x)+y, -x, z; 
%           y, x, (-z)+1/2; x-y, -y, (-z)+1/2; -x, (-x)+y, (-z)+1/2;
%           -x, -y, -z; y, (-x)+y, -z; x-y, x, -z; 
%           -y, -x, z+1/2; (-x)+y, y, z+1/2; x, x-y, z+1/2];
%    names_36f = {'36f', '36f', '36f', '36f', '36f', '36f', '36f', '36f', ...
%        '36f', '36f', '36f', '36f'};
  %All the positions in this material     
    pos_sym = [pos_6b(:,:); pos_18e(:,:)]; %pos_36f(:,:)];
    names_sym = {names_6b{:}, names_18e{:}}; %names_36f{:}};   
  % Construct all:  
   [names, pos] = construct(sym', pos_sym', names_sym, supercell);
end
  
%% Li3PS4-beta
function [names, pos] = li3ps4_beta(supercell)
%The symmetries applicable to each site:
    sym = [0 0 0];
    % 8d positions
    % Homma et al. (2011): 
    %x = 0.3562; y = 0.013; z = 0.439; 
    % Mercier et al. (1982): 
    %x = 0.332; y = 0.0335; z = 0.3865;
    % Yan Chen et al. (2015), Appl. Phys. Lett.
    x = 0.341; y = 0.031; z = 0.384; %(413 K)
    
    pos_8d = [ x y z; (-x)+0.5 -y z+0.25; -x y+0.5 -z; x+0.5 (-y)+0.5 (-z)+0.25; 
        -x -y -z; x+0.5 y (-z)+0.25; x (-y)+0.5 z; (-x)+0.5 y+0.5 z+0.25 ];
    names_8d = {'8d', '8d', '8d', '8d', '8d', '8d', '8d', '8d'};
    
    %4b positions
    pos_4b = [ 0 0 0.5; 0.5 0 0; 0 0.5 0.5; 0.5 0.5 0 ];
    names_4b = {'4b' '4b' '4b' '4b'};
    
    %4c positions
    % Homma et al. (2011): 
    %x = 0.999; z = 0.889; 
    % Mercier et al. (1982): 
    %x = -0.074; z = -0.306; 
    % Yan Chen et al. (2015), Appl. Phys. Lett.
    x = -0.017;  z = -0.263; %annealed at 413 K
    %x = -0.056; z = -0.361; %annealed at 673 K
    pos_4c = [x 0.25 z; (-x)+0.5 0.75 z+0.25; -x 0.75 -z; x+0.5 0.25 (-z)+0.25 ];
    names_4c = {'4c' '4c' '4c' '4c'};

 % All the positions in this material     
    pos_sym = [pos_8d(:,:); pos_4b(:,:); pos_4c(:,:)];
    names_sym = {names_8d{:}, names_4b{:}, names_4c{:}};   
  % Construct all:  
    [names, pos] = construct(sym', pos_sym', names_sym, supercell);
   
   disp('Giving different names to detect intersheet jumps')
   for i = 1:size(pos,2)
        if pos(1,i) > 0.25 && pos(1,i) < 0.75
            names{i} = [names{i}, 'sheet1'];
        else
            names{i} = [names{i}, 'sheet2'];
        end
   end
   
end

function [names, pos] = construct(sym, pos_sym, names_sym, supercell)
% Constructs all the (fractional) positions based on the symmetries and coordinates,
% and combines it with the given names.

%% Check whether names_sym and pos_sym have the same length:
    if size(pos_sym,2) ~= size(names_sym)
        disp('ERROR! Positions and names given not the same size!')
        size(pos_sym,2)
        size(names_sym,2)
        %names_sym
    end
    
% The number of final positions:
    nr_sym = size(sym,2);
    pos_given = size(pos_sym,2);
    nr_cells = prod(supercell);
    pos = zeros(3, nr_cells*nr_sym*pos_given);
    names = cell(nr_cells*nr_sym*pos_given,1);
    counter = 0;
    
%% Start constructing the unit cell based on the given symmetry operations and coordinates:
    for i = 1:pos_given
        for j = 1:nr_sym %the symmetries           
            counter = counter + 1;
            names{counter} = names_sym{i};
            for k = 1:3 %the a,b,c coor  
                pos(k, counter) = pos_sym(k,i) + sym(k,j);
                if pos(k, counter) >= 1 % Periodic Boundary Conditions:
                    pos(k, counter) = pos(k, counter) -1;
                elseif pos(k, counter) < 0
                    pos(k, counter) = pos(k, counter) + 1;
                end
                % Take coordinate shift of fractional coordinates due to supercell into account
                if supercell(k) > 1
                    pos(k, counter) = pos(k, counter)/supercell(k);
                end
            end
        end
    end
    
%% Construct the supercells:
    if nr_cells > 1
        for a = 1:supercell(1)
            for b = 1:supercell(2)
                for c = 1:supercell(3)                
                    if a ~= 1 || b ~= 1 || c ~= 1 
                    % The [1 1 1] cell has already been constructed, now do only do the rest   
                        for i = 1:nr_sym*pos_given
                            counter = counter + 1;
                            names{counter} = names{i};
                            pos(1,counter) = pos(1,i) + (a-1.0)/supercell(1);
                            pos(2,counter) = pos(2,i) + (b-1.0)/supercell(2);
                            pos(3,counter) = pos(3,i) + (c-1.0)/supercell(3);
                        end

                   end
                end %c 
           end
       end
    end
   % Test: scatter3(pos(1,:), pos(2,:), pos(3,:))
    
end