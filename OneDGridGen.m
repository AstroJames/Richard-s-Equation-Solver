function [Del,z,hInitial] = OneDGridGen(delt,n,L,viz,type)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OneDGridGen produces a one dimensional grid that can either be uniform or
% take on geometric properties depending on the user input. The function is
% also able to visualise the grid that has been produced. The function
% outputs all of the changes in time, cells and nodes as a cell array.

% Inputs:
% X     = the size of the grid. 
% DELT  = the time step.
% N     = the number of points on the grid.
% VIZ   = turns the grid visualisation on.
% TYPE  = the type of grid: uniform, left geometric, right geometric,
%           bilateral geometric
%
%
% Outputs:
% DEL       = a cell array of all of the steps in space and time. 
% Delx      = differences in cells
% delx      = differences in nodes
% delt      = differences in time
% hInitial  = the initial guess of the pressure head
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate the difference between nodes and cells


switch type
    % Build a uniform grid.
    case 'uniform'
        disp('You have selected uniform for the grid distribution.')
        z    = 0:L/(n-1):L;
    % Create a geometrically distributed left grid. 
    case 'left geometric'
        % Initialisations
        disp('You have selected left geometric for the grid distribution.')
        r               = 1.05;
        geometric_x     = L*(1-r) / (1-r^n); 
        z               = zeros(1,n);
        
        % Create the grid
        for i=1:n-1
            z(i+1)   = z(i) + r^(i-1)*geometric_x;
        end
    % Create a geometrically distributed right grid.    
    case 'right geometric'
        % Initialisations
        disp('You have selected right geometric for the grid distribution.')
        r               = 1.05;
        geometric_x     = L*(1-r) / (1-r^n); 
        z               = zeros(1,n);  
        
        % Create the grid
        for i=1:n-1
            z(i+1)   = z(i) + r^(n-i)*geometric_x;
        end
    % Create a geometrically distributed right grid.
    case 'bilateral geometric'
        % Initialisations
        disp('You have selected bilateral geometric for the grid distribution.')
        r               = 1.05;
        z               = zeros(1,n);
        k1              = length(1:round(n/2));
        k2              = length((n-round(n/2)):n-1);
        
        geometric_x     = round(L/2)*(1-r) / (1-r^k1);
            
        % Start on one side and build a geometric grid
        for i=1:round(n/2)
            z(i+1)   = z(i) + r^(i-1)*geometric_x;
        end
        
        % Start from the other side and build a geometric grid
        geometric_x     = (L-round(L/2))*(1-r) / (1-r^k2); 
        
        for i=(n-round(n/2)):n-1
            z(i+1)   = z(i) + r^(n-i)*geometric_x;
        end 
    otherwise 
        error('Error. Please choose the grid distribution. Your options are `uniform`, `left geometric`, `right geometric` or `bilaterally geometric`.')
end

delx = diff(z);             % difference between nodes
Delx = zeros(1,n);          % difference between cells

% Create the distance between the cells.
Delx(1) = (z(2)-z(1))/2;
Delx(2:n-1)= (z(3:n)-z(1:n-2))/2;
Delx(n) = (z(n)-z(n-1))/2;



%% Visualisation of the Nodes

% Build a rudimentary visualisation of the placement of the nodes and cells
% along the grid.
switch viz
    case 'on'
        % Find the placement of the cells
        cells = zeros(1,length(Delx));
        cells(1) = Delx(1);
        for i = 2:length(Delx)
            if sum(i == 2:length(Delx)-1) == 1
                cells(i) = sum(Delx(1:i));
            else
            cells(i) = sum(Delx(1:i-1))+2*Delx(length(Delx)); 
            end
        end
        % plot nodes
        plot(z,0,'.','color','black','MarkerSize',20)
        hold on
        % plot cells
        plot(cells,0,'o','color','red','MarkerSize',5)
        title('Nodes indicated in black, Cell boundaries indicated in red')
    case 'off'
        disp('You have selected for no visualisation of the grid.')
    otherwise
        error('Error. Please choose either grid visualisation on, `on`, or off `off`')
        
end


% All steps in space and time.
Del         = {Delx',delx',delt};
z           = z'; 

% The initial guess of the pressure head.
hInitial    = zeros(n,1);
for i = 1:n
    hInitial(i) = -z(i) + 70; 
end

end % end function