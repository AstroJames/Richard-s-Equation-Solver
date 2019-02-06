function [x,z,del,m,n] = twod_gridgenerator(m, n, type,viz,r,s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% twod_gridgenerator produce a 2D grid containing the discretised nodes of
% the 2D homogenisation problem, as well as any relevant information
% (distance between nodes, control volumes, etc.).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% m         =   the number of nodes in the x-direction.
% n         =   the number of nodes in the z-direction.
% type      =   the type of grid spacing distribution, where:
%               - 'uniform' use uniformly space the grid such that the unit
%               cell can be evenly subdivided into 4x4 cells,
%               - 'random' generates random grid spacing through uniform
%               probability distribution, where [0,0.25,0.5,0.75,1] are
%               always a coordinate of the nodes in x and z directions,
%               - 'bilateral' generates a grid using bilateral geometric
%               spacing in each 0.25 subgrid, with the spacing ratio of
%               s.
% viz       =   a string containing 'on' or 'off', which determines the
%               inclusion of visualisation of the 2D Grid Generator.
% Optional Inputs:
% r         =   the type of the material of interest, where
%               - 'a' is the Alluvium Aquifer,
%               - 'c' is the Walloon Coal Measure,
%               - 't' is the benchmark test case.
%               If specified, visualisation will be limited to the type of
%               material specified.
%               If 'a' is selected when type is 'bilateral', the grid
%               deploy special bilateral distribution specifically designed
%               for Alluvium Aquifer unit cell.
%               If not specified, the value is default to 'c', but
%               visualisation will display both Alluvium Aquifer and
%               Walloon Coal unit cells.
% s         =   the grid spacing ratio in x and z direction. If not
%               specified, this input is default to 2.
% Outputs:
% x         =   a vector containing the coordinates of each nodes in the
%               x-direction.
% z         =   a vector containing the coordinates of each nodes in the
%               z-directions
% del       =   a cell containing 2 vectors in the following order:
%               {delx, delz}, where
%               - delx is a vector containing the x-distance between nodes,
%               - delz is a vector containing the z-distance between nodes.
% m         =   the modified number of nodes as considered in the
%               x-direction.
% n         =   the modified number of nodes as considered in the
%               z-direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Due to the alluvium bilateral grid being modified for Group 1
% unit cell geometry specifically, having 'bilateral' grid type and 'a'
% materal type will most-likely be incompatible with other unit cells
% geometry, and may break the code.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check for optional inputs
ncheck = nargin;
if nargin < 6
    if nargin == 5 && ischar(r) % r input only
        s = 2;
        ncheck = nargin + 1;
    elseif nargin == 5 && ~isstring(r) % s input only
        s = r;
        r = 'c';
    else % No optional inputs
        s = 2;
        r = 'c';
    end
end
% Check for other inputs
if  isreal(m) == 0 || isreal(n) == 0 || mod(m,1) ~= 0 || mod(n,1) ~= 0
    error('The number of nodes must be real positive integers')
end
if m < 5 || n < 5
    error('The number of nodes in x and z must be 5 or more.')
end

%% Grid generator

switch type
    case 'uniform'
        % Check and modify number of nodes
        m_new = 4*round((m-1)/4)+1;
        if m_new ~= m
            m = m_new;
            disp('The value of m has been changed due to compatability issue.')
        end
        
        n_new = 4*round((n-1)/4)+1;
        if n_new ~= n
            n = n_new;
            disp('The value of n has been changed due to compatability issue.')
        end
        
        % Uniform grid generation
        x = 0:1/(m_new-1):1;
        z = 0:1/(n_new-1):1;
        
    case 'random'
        % Defining constraints
        v = [0,0.25,0.5,0.75,1];
        
        % Random grid generation
        x = sort([rand(1,m-5),v]);
        z = sort([rand(1,n-5),v]);
        
        % Safeguard function
        x = repeat_remover(x);
        z = repeat_remover(z);
    case 'bilateral'
        % Check and modify number of nodes
        m_local = round((m-1)/4);
        m_new = 4*m_local+1;
        if m_new ~= m
            m = m_new;
            disp('The value of m has been changed due to compatability issue.')
        end
        
        n_local = round((n-1)/4);
        n_new = 4*n_local+1;
        if n_new ~= n
            n = n_new;
            disp('The value of n has been changed due to compatability issue.')
        end
        
        % Bilateral geometric grid generator
        if r == 'a'
            x = modified_bilateral_generator(m_local,s);
        else
            x = bilateral_generator(m_local,s);
        end
        z = bilateral_generator(n_local,s);
    otherwise
        error('Error: Invalid grid distribution "type".')
end

%% Calculation of Del

delx = diff(x);                 % Distance between nodes in x - direction
delz = diff(z);                 % Distance between nodes in z - direction
del = {delx, delz} ;            % Grouping output

% Visualisation
switch viz
    case 'on'
        disp('You have select visualisation on for the 2D grid generator.')
        
        cells_x = (x(2:m)+x(1:m-1))/2;  % coordinates of the control volume face
        cells_z = (z(2:n)+z(1:n-1))/2;  % coordinates of the control volume face
        cells = {cells_x, cells_z};     % Grouping cells
        
        if ncheck == 6 % If material type has been specified
            gridandk_visualisation(x,z,cells,r)
        else % If material type has not been specified
            gridandk_visualisation(x,z,cells,'a')
            gridandk_visualisation(x,z,cells,'c')
        end
        
    case 'off'
        disp('You have select visualisation off for the 2D grid generator.')
    otherwise
        error('Please select "on" or "off" for visualisation.')
end

end

function y = repeat_remover(x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In a random grid case, repeat_remover is a safeguard function designed to
% remove nodes that may overlap each other while also ensuring that the
% number of nodes doesn't change.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% x     =   the original coordinates of the grid.
% Output:
% y     =   the modified coordinates of the grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialising variables
y = x;
n = length(y);
a = 1; % Initialised so the safeguard run at least once

% Safeguard function
while a < n % Check that the length of the vector doesn't change
    y = sort([y,rand(1,n-length(y))]); % Compensate for length change if any
    i = 1;
    while i < length(y)
        if y(i) == y(i+1) % Check for repeated values
            y = [y(1:i),y(i+2:end)]; % Remove the said repeated value
        else
            i = i+1;
        end
    end
    a = length(y);
end

end

function gridandk_visualisation(x,z,cells,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% gridandk_visualisation shows the visualisation of the discretised nodes
% on the unit cells. It also compare the value of K at each node against
% the material types of the unit cell where the node is located.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% x         =   a vector containing the coordinates of each nodes in the
%               x-direction.
% z         =   a vector containing the coordinates of each nodes in the
%               z-directions
% cells     =   a cell containing 2 vectors in the following order:
%               {cells_x, cells_z}, where
%               - cells_x is a vector containing the x-coordinates of the
%               control volume boundary (excluding the unit cell boundary).
%               - cells_z is a vector containing the z-coordinates of the
%               control volume boundary (excluding the unit cell boundary).
% r         =   the type of the material of interest, where
%               - 'a' is the Alluvium Aquifer
%               - 'c' is the Walloon Coal Measure
%               - 't' is the benchmark test case.
% Output:
% plot a visualisation of the discretisation of the unit cell, as well as
% check the value of K on each note against the material type.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise variables
cells_x = cells{1};
cells_z = cells{2};
soilA = 1/255*[204 229 255];

% Initialise figure
figure
hold on
% Start visualisation
switch r
    case 'a' % Alluvium Aquifer
        % Set up material A
        rectangle('position', [0.25, 0.25, 0.75, 0.25], 'facecolor', soilA)
        rectangle('position', [0, 0.75, 0.75, 0.25], 'facecolor', soilA)
        % K visualisation
        for i = 1:length(x)
            for j = 1:length(z)
                if ((z(j)==0.25||z(j) == 0.5)&&x(i)>=0.25)||((x(i)==0.25||mod(x(i),1) ==0)&&(z(j)>=0.25 && z(j)<=0.5))
                    plot(x(i),z(j), '.','color', 'k', 'MarkerSize', 15) % Lower blob boarder
                elseif ((x(i)==0.75||mod(x(i),1)==0)&&(z(j)>=0.75))||((z(j)==0.75||mod(z(j),1)==0)&&(x(i)<=0.75))
                    plot(x(i),z(j), '.','color', 'k', 'MarkerSize', 15) % Upper blob boarder
                elseif x(i) == 1 && z(j) == 0
                    plot(x(i),z(j), '.','color', 'k', 'MarkerSize', 15) % Blindspot boarder
                elseif K(x(i),z(j),r) == 6.07
                    plot(x(i),z(j), '.','color', 'b', 'MarkerSize', 15) % Material A
                else
                    plot(x(i),z(j), '.','color', 'r', 'MarkerSize', 15) % Material B
                end
            end
        end
    case 'c' % Walloon Coal
        % Set up material A
        for i = 0:1
            for j = 0:1
                rectangle('position', [0.5*i, 0.5*j, 0.25, 0.25], 'facecolor', soilA)
                rectangle('position', [0.5*i+0.25, 0.5*j+0.25, 0.25, 0.25], 'facecolor', soilA)
            end
        end
        % K visualisation
        for i = 1:length(x)
            for j = 1:length(z)
                if mod(x(i),0.25) == 0 || mod(z(j),0.25) == 0
                    plot(x(i),z(j), '.','color', 'k', 'MarkerSize', 15) % Boarders
                elseif K(x(i),z(j),r) == 2.32*10^-1
                    plot(x(i),z(j), '.','color', 'b', 'MarkerSize', 15) % Material A
                else
                    plot(x(i),z(j), '.','color', 'r', 'MarkerSize', 15) % Material B
                end
            end
        end
    case 't' % Benchmark test case. 
        % Set up material A
        rectangle('position', [0, 0, 1, 0.5], 'facecolor', soilA)
        % K visualisation
        for i = 1:length(x)
            for j = 1:length(z)
                if mod(z(j),0.5) == 0
                    plot(x(i),z(j), '.','color', 'k', 'MarkerSize', 15) % Boarders
                elseif K(x(i),z(j),r) == 16
                    plot(x(i),z(j), '.','color', 'b', 'MarkerSize', 15) % Material A
                else
                    plot(x(i),z(j), '.','color', 'r', 'MarkerSize', 15) % Material B
                end
            end
        end
    otherwise
        error('Unknown input material type.')
end

% Grid visualisation
line([0;1],[cells_z;cells_z],'color','r') % Draw horizontal lines
line([cells_x;cells_x],[0;1],'color','r') % Draw vertical lines
line([0,0,1,1,0],[0,1,1,0,0],'color','b') % Draw boarders
axis equal

end

function x = bilateral_generator(m,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bilateral_generator creates a vector containing the coordinates with
% bilateral geometric spacing for each sub-grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% m         =   the number of nodes as considered in the x-direction.
% r         =   the spacing ratio for the geometric spacing.
% Output:
% x         =   a vector containing the coordinates of each nodes in the
%               x-direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The direction defined above is only valid within this function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialised grid
x = zeros(1,4*m+1);
x([m+1,2*m+1,3*m+1,4*m+1]) = [0.25,0.5,0.75,1]; % Forced of subgrids boundary

% Spacing of the sub-grid
if mod(m,2) == 0 % Even case
    k = 0:m/2-2;
    del_base = (1-r)/(8*(1-r^(m/2))); % Calculate base distance
    del_local = del_base*r.^k; % Calculate spacing
    x([m/2+1,3*m/2+1,5*m/2+1,7*m/2+1]) = [1/8,3/8,5/8,7/8]; % Forced midpoint of subgrids
else % Odd case
    k = 0:(m+1)/2-2;
    del_base = (1-r)/(4*(2-r^((m+1)/2-1)*(1+r))); % Calculate base distance
    del_local = del_base*r.^k; % Calculate spacing
end
% Grid construction
for j = 1:length(del_local)
    % Bottom half of subgrids
    x([j+1,m+j+1,2*m+j+1,3*m+j+1]) = x([j,m+j,2*m+j,3*m+j]) + del_local(j);
    % Top half of subgrids
    x([m+1-j,2*m+1-j,3*m+1-j,4*m+1-j]) = x([m+2-j,2*m+2-j,3*m+2-j,4*m+2-j]) - del_local(j);
end
end

function x = modified_bilateral_generator(m,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% modified_bilateral_generator creates a vector containing the x-coordinate
% with a special bilateral geometric spacing designed specifically for the
% Alluvium Aquifer's unit cell.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% m         =   the number of nodes as considered in the x-direction.
% r         =   the spacing ratio for the geometric spacing.
% Output:
% x         =   a vector containing the coordinates of each nodes in the
%               x-direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialised grid
x = zeros(1,4*m+1);
x([m+1,3*m+1,4*m+1]) = [1/4,3/4,1]; % Forced of subgrids boundary

% Spacing of the corner sub-grid
if mod(m,2) == 0 % Even case
    k = 0:m/2-2;
    del_base = (1-r)/(8*(1-r^(m/2))); % Calculate base distance
    del_local = del_base*r.^k; % Calculate spacing
    x([m/2+1,7*m/2+1]) = [1/8,7/8]; % Forced midpoint of subgrids
else % Odd case
    k = 0:(m+1)/2-2;
    del_base = (1-r)/(4*(2-r^((m+1)/2-1)*(1+r))); % Calculate base distance
    del_local = del_base*r.^k; % Calculate spacing
end

% Spacing of centre sub-grid
k_centre = 0:m-2;
del_base_centre = (1-r)/(4*(1-r^(m))); % Calculate base distance
del_local_centre = del_base_centre*r.^k_centre; % Calculate spacing
x(2*m+1) = 1/2; % Forced midpoint of subgrids

% Grid construction - corner sub-grid
for j = 1:length(del_local)
    % Bottom half of subgrids
    x([j+1,3*m+j+1]) = x([j,3*m+j]) + del_local(j);
    % Top half of subgrids
    x([m+1-j,4*m+1-j]) = x([m+2-j,4*m+2-j]) - del_local(j);
end
% Grid construction - centre sub-grid
for j = 1:length(del_local_centre)
    % Bottom half of subgrids
    x(m+j+1) = x(m+j) + del_local_centre(j);
    % Top half of subgrids
    x(3*m+1-j) = x(3*m+2-j) - del_local_centre(j);
end
end