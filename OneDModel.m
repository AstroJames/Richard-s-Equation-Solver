function [hStorage,runTime,timeSteps,FEvals,epsis] = OneDModel(tInterval,z,Del,Kzz,Kxx,hInitial,CSGPeriods,droughtPeriods,CSGDistance,rainSetting,solveSetting,VariableCSG,m,sigma,theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is the main driving function for the 1D Model
%
% INPUTS:
% tInterval         = is the entire time period for the model,
%                       tInterval(1) is start,
%                       tInterval(2) is end.
% z                 = the vector that contains all of the height of all the nodes.
% Del               = the cell that contains all the differences,
%                       Del{1} is volume of each cell,
%                       Del{2} is distance between each node,
%                       Del{3} is the initial time step
% Kzz               = the vector that conatians the Kzz values for the different
%                     sections in the model,
%                       Kzz(1) is aqufier,
%                       Kzz(2) is coal measure.
% Kxx               = the vector that conatians the Kxx values for the different
%                     sections in the model,
%                       Kxx(1) is aqufier,
%                       Kxx(2) is coal measure.
% hInitial          = the vector containing initial pressure for all nodes.
% CSGPeriods        = the matrix that contains all the periods of when the CSG
%                     plant is on,
%                       1st column is start times,
%                       2nd column is end times.
% droughtPeriods    = the matrix that contains all the periods of when the
%                     model is in drought,
%                       1st column is start times,
%                       2nd column is end times.
% CSGDistance       = the distance of CSG Plant.
% rainSetting       = the setting of the rain function, takes:
%                       'High Frequency', 'Low Frequency', 
%                       'Low Frequency Stochastic', and 'Simple'
% solveSetting      = string for how the solver will find the Newton
%                     step, takes
%                        'TDMASolve' and 'GMRES'
% VariableCSG       = the string that is the condition of whether the CSG plant
%                     is varrying, takes:
%                        'on' and 'off'
% m                 = Newton settings. 0 is Chord, 1 is Full and, 1< is Shamanski
% sigma             = the control on what aproximation will be used on the velocity
%                     terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%                     upwinding
% theta             = setting for the interation,
%                        theta(1) is for flux,
%                        theta(2) is for source.
% OUTPUTS:
% hStorage          = the matrix containing the presures for all node over
%                     all timesteps.
% runTime           = the total running time for the model.
% timeSteps         = the vector containing the amount of time incroment 
%                     for each time step.
% FEvals            = the vector containing the number of function
%                     evaulations each time step
% epsis             = the vector containing the water buget error each time
%                     step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Input Checks
if min(min(CSGPeriods)) < tInterval(1) || min(min(droughtPeriods)) < tInterval(1)
    error('One of the inputed periods for CSG Extraction or Droughts starts before the time interval')
end

if max(max(CSGPeriods)) > tInterval(2) || max(max(droughtPeriods)) > tInterval(2)
    error('One of the inputed periods for CSG Extraction or Droughts ends after the time interval')
end

%% Settings for the Model
tol               = [0,10e-8];        % relative, absolute tolerance
lineSearchSetting = 'two-point';


%% Variable Initialization
h               = [hInitial,hInitial];  % the initial guess for h 
zAt51           = FindZAt(z,51);        % the place in the grid which corresponds z = 51
R               = [0;0];                % initialise the weather function output
t               = tInterval(1);         % intialise the time
i               = 1;                    % intialise the total number of iterations through the time stepping
conInRowCounter = 1;                    % intialise the number of convergences in the row to update the time step
CSGCounter      = 0;                    % initialise the counter for the number of convergences in a row.
droughtCounter  = 0;                    % the drought state is set to off
FEvals          = zeros(1,length(tInterval(1):Del{3}:tInterval(2))); % an initialisation of the number of F(x)


%% Control Constants
jMax           = 10;    % the maximum for trying the Newton method
deltMax        = 1;     % the maxmum time step = 1 day
conInRowNumber = 5;     % Number of Sucessful Iteration before increasing time step

%% Storage Setup
hStorage       = zeros(size(h,1),length(tInterval(1):Del{3}:tInterval(2))); % intialising the pressure head matrix
hStorage(:,1)  = h(:,2);                                                    % storing the first guess
timeSteps      = zeros(1,length(tInterval(1):Del{3}:tInterval(2)));         % initialising the vector of time steps the solver takes
epsis          = zeros(1,length(tInterval(1):Del{3}:tInterval(2)));         % initialising the vector associated with the water budget

%% Running the Model
tic
while t < tInterval(2)
    
    j       = 0; 
    converg = false;
    fprintf('The current iteration is %d. \n',i)
    
    % Newton Solver
    while converg == false && j < jMax
        dt          = [t; t+Del{3}];
        
        % Drought Handling
        if droughtPeriods(1,1) ~= droughtPeriods(1,2)
            if droughtCounter < size(droughtPeriods,1)
                droughtState = StateFinder(dt,droughtPeriods(droughtCounter+1,:));
            end
        else
            droughtState = 'FullyOut';
        end
        
        
        % CSG Handling
        if CSGPeriods(1,1) ~= CSGPeriods(1,2)
            if CSGCounter < size(CSGPeriods,1)
                CSGState = StateFinder(dt,CSGPeriods(CSGCounter+1,:));
            end
        else
            CSGState = 'FullyOut';
        end
        
        %Rainfall
        if t == 0 
            [R(1),R(2)] = SourceFunction(rainSetting,t,Del,droughtState);
            Rnull       = R(2);
        else
            [~,R(2)]    = SourceFunction(rainSetting,t,Del,droughtState);
            R(1)        = Rnull;
            Rnull       = R(2);
        end
        
        % where all the work happens
        try
            [h,epsis(i),FEvals(i)]       = NewtonMethod(h,dt,z,Del,R,Kzz,Kxx,zAt51,CSGState,droughtState,CSGDistance,theta,sigma,tol,m,lineSearchSetting,solveSetting,VariableCSG);
            converg             = true;
        catch ErrorMes
            switch ErrorMes.identifier
                case 'LineSearch:NoConverg'
                    Del{3}    = Del{3}/2;
                    conInRowCounter    = 1;
                    fprintf('The time step has been shrunk to %d. \n',Del{3})                    
                otherwise
                    rethrow(ErrorMes)
            end
        end
        
        j   = j + 1;
    end
    
    % Convergence test
    if j == jMax
        error('Timestep halfing has reached the max iterations')
    elseif j == 1
        conInRowCounter = conInRowCounter + 1;
    end
    
    % Store Solution
    h(:,2)              = h(:,1);
    hStorage(:,i+1)     = h(:,1); 
    
    % Make time step larger after conInRowNumber good iterations have been made 
    if conInRowCounter == conInRowNumber
        if Del{3} < deltMax
            Del{3}              = Del{3}*1.1; 
            conInRowCounter     = 1; 
            fprintf('The time step has been grown to %d. \n',Del{3})
        else
            Del{3}    = deltMax;
        end
    end
    
    timeSteps(i)     = Del{3};
    
    
    
    i               = i + 1; 
    t               = t + Del{3};
    
    %Updating counters for how many droughts and CSG extraction periods
    %have occured
    if droughtCounter < size(droughtPeriods,1) &&t > droughtPeriods(droughtCounter+1,2)
        droughtCounter = droughtCounter+1;
    end
    
    if CSGCounter < size(CSGPeriods,1) && t > CSGPeriods(CSGCounter+1,2)
        CSGCounter = CSGCounter+1;
    end
end
runTime = toc;
%Clearing all the outputs
timeSteps(~any(timeSteps,1)) = [];
epsis(~any(epsis,1)) = [];
FEvals(~any(FEvals,1)) = [];
hStorage( :, ~any(hStorage,1) ) = [];

end


%% Extra Functions

function zAtNum = FindZAt(z,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the node closest to a particular height.
%
% INPUTS:
% z         = the height of all nodes.
% num       = a particular height.
% OUTPUTS:
% zAtNum    = a vecto that contains all the infomation of the closed node at
%             num,
%               zAtNum(1) is the index of where the node is,
%               zAtNum(2) is the actual node height.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zAtNum      = zeros(2,1);
zAtNum(1)   = 1;
check       = inf;

while abs(num - z(zAtNum(1))) < check
    check       = num - z(zAtNum(1)); %Reseting the closet distance
    zAtNum(1)   = zAtNum(1)+1; % next node to test
end

zAtNum(2)   = z(zAtNum(1)); 
end

function state = StateFinder(tInterval,timePeriod)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function the the particular state the model is going into. 
%
% INPUTS:
% tInterval     = a vector containing the current times in the model.
%                   tInterval(1) is the forward part,
%                   tInterval(2) is the backward part,
% timePeriod    = the next time period in the model.
%                   timePeriod(1) is the start of the period,
%                   timePeriod(2) is the end of the period,
% OUTPUTS:
% state         = the string that tells what state the model is in, throws:
%                   'FullyIn', 'GoingOut', 'FullyOut' and 'GoingIn'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if tInterval(1) >= timePeriod(1)
    if tInterval(2) <= timePeriod(2)
    	state = 'FullyIn';
    else
    	state = 'GoingOut';
	end
else
	if tInterval(2) < timePeriod(1)
    	state = 'FullyOut';
    else
    	state = 'GoingIn';
	end
end

end