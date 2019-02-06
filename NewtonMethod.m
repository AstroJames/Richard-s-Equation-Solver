function [h,epsi,FCounter] = NewtonMethod(h,dt,z,Del,R,Kzz,Kxx,zAt51,CSGState,droughtState,CSGDistance,theta,sigma,tol,m,lineSearchSetting,solveSetting,VariableCSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function that handles the Newtons method for the 1D model.
%
% INPUTS:
% h                    = the vector of all the pressure heads at each node at each 
%                        time step.
% dt                   = is the time vector,
%                           dt(1) is forward time step,
%                           dt(2) is backward time step.
% z                    = the vector that contains all of the height of all the nodes.
% Del                  = the cell that contains all the differences,
%                           Del{1} is volume of each cell,
%                           Del{2} is distance between each node,
%                           Del{3} is the time difference.
% R                    = the rainfall on the suface (m)
% Kzz                  = the vector that conatians the Kzz values for the different
%                        sections in the model,
%                           Kzz(1) is aqufier,
%                           Kzz(2) is coal measure.
% Kxx                  = the vector that conatians the Kxx values for the different
%                        sections in the model,
%                           Kxx(1) is aqufier,
%                           Kxx(2) is coal measure.
% zAt51                = vector that contains infomation above the node closest to the
%                        height z = 51,
%                           zAt51(1) is actual height at the node,
%                           zAt51(2) is the index of where the node is.
% CSGstate             = the state of CSG extraction, takes:
%                           'FullyIn', 'GoingOut', 'FullyOut', and 'GoingIn'
% droughtState         = the state of drought in the model, takes:
%                           'FullyIn', 'GoingOut', 'FullyOut', and 'GoingIn'
% CSGDistance          = the distance of CSG Plant.
% theta                = setting for the interation,
%                           theta(1) is for flux,
%                           theta(2) is for source.
% sigma                = the control on what aproximation will be used on the velocity
%                        terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%                        upwinding
% tol                  = tollerance of error in the solver. 
%                           tol(1) relative tolerance
%                           tol(2) absolute tolerance
% m                    = Newton settings. 0 is Chord, 1 is Full and, 1< is Shamanski
% lineSearchingSetting = string for what type of linesearching, takes:
%                           'simple' and 'two-point'
% solveSetting         = string for how the solver will find the Newton
%                        step, takes
%                           'TDMASolve' and 'GMRES'
% VariableCSG          = the string that is the condition of whether the CSG plant
%                        is varrying, takes:
%                           'on' and 'off'
% OUTPUTS:
% h                    = the vector of all the pressure heads at each node at each 
%                        time step.
% espi                 = error in the water bugdet.
% FCounter             = number of times F is computed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Error checks 

%% Solves the CSG 1D system using inexact Newton's Method
FCounter = 0;        
i = 0;
maxIt = 5;
CSGOn           = DeterminState(CSGState);
droughtOn       = DeterminState(droughtState);
[Fh0,FCounter]  = Ffunc(h,dt,z,zAt51,Del,R,Kzz,Kxx,CSGOn,droughtOn,CSGDistance,theta,sigma,FCounter,VariableCSG);
%FhOld           = Fh0;
Fh              = Fh0;
x0              = zeros(length(h),1);       % initial guess for gmres implementation.
maxiters        = 500;                      % maximum iterations for gmres implementation.


while norm(Fh,2) > tol(1)*norm(Fh0,2) + tol(2) && i < maxIt
    switch solveSetting
        case 'TDMASolve'
            %Calculate Jacobian
            if mod(i,m) == 0 %|| norm(Fh)/norm(FhOld) > 0.95
                [J,FCounter] = JacobianTri(h,dt,z,zAt51,Del,R,Kzz,Kxx,CSGOn,droughtOn,CSGDistance,theta,sigma,Fh,FCounter,VariableCSG);
            end  
            
            delh = - TDMA_solve(J,Fh); %Calculate Newton Step
        case 'GMRES'
            %Calculate Jacobian
            [J,FCounter]     = JacobianTri(h,dt,z,zAt51,Del,R,Kzz,Kxx,CSGOn,droughtOn,CSGDistance,theta,sigma,Fh,FCounter,VariableCSG);
            if mod(i,10) == 0
                JOld = J; %Store as Old if hit 10 iteration
            end
            delh = - generalmres_oned(x0, J, Fh, tol(2), maxiters, JOld); %Calculate Newton Step
        otherwise
            error('You need to select either GMRES or TDMASolve to select the Newton solver.')
    end
    
    %Find h using line searching
    [h,Fh,FCounter] = LineSearch(h,dt,Fh,z,zAt51,Del,R,Kzz,Kxx,CSGOn,droughtOn,CSGDistance,theta,sigma,delh,lineSearchSetting,FCounter,VariableCSG);
    %FhOld  = Fh;
    i = i + 1; %increment counter
end

epsi = EpsiCal(h,dt,z,zAt51,Del,Kxx,CSGOn,droughtOn,R,CSGDistance,VariableCSG);

end


function [J,FCounter] = JacobianTri(x,dt,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,Fx,FCounter,VariableCSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JacobianTri finds the Tridigaonal Jacobian for fFunc at x
%
% INPUTS:
% x           = the vector of all the pressure heads at each node at each 
%               time step.
% dt          = is the time vector,
%                   dt(1) is forward time step,
%                   dt(2) is backward time step.
% z           = the vector that contains all of the height of all the nodes.
% zAt51       = vector that contains infomation above the node closest to the
%               height z = 51,
%                   zAt51(1) is actual height at the node,
%                   zAt51(2) is the index of where the node is.
% Del         = the cell that contains all the differences,
%                   Del{1} is volume of each cell,
%                   Del{2} is distance between each node,
%                   del{3} is the time difference.
% R           = the rainfall on the suface (m)
% Kzz         = the vector that conatians the Kzz values for the different
%               sections in the model,
%                   Kzz(1) is aqufier,
%                   Kzz(2) is coal measure.
% Kxx         = the vector that conatians the Kxx values for the different
%               sections in the model,
%                   Kxx(1) is aqufier,
%                   Kxx(2) is coal measure.
% CSGon       = the state of CSG extraction, 1 is on and 0 is off.
% droughtOn   = the state of drought in the model, 1 is in drough and 0 is
%               out of drought.
% CSGDistance = the distance of CSG Plant.
% theta       = setting for the interation,
%                   theta(1) is for flux,
%                   theta(2) is for source.
% sigma       = the control on what aproximation will be used on the velocity
%               terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%               upwinding
% Fx          = the function evaluated at x.
% FCounter    = counts of how many time F is computed.
% VariableCSG = the string that is the condition of whether the CSG plant
%               is varrying, takes:
%                   'on' and 'off'
% OUTPUTS:
% J           = The Jacobian for the new system.
% FCounter    = counts of how many time F is computed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = length(x);
J{1} = zeros(n-1,1);
J{2} = zeros(n,1);
J{3} = zeros(n-1,1);

%Set up s vectors
s = zeros(n,3);
for i = 0:n-1
    s(i+1,mod(i,3)+1) = 1;    
end

%Compute components for J
for i = 1:3
    if norm(x,2) ~= 0
        h = sqrt(eps)*norm(x,2)/norm(s(:,i),2);
    else
        h = sqrt(eps)/norm(s(:,i),2);
    end
    s(:,i) = ( Ffunc([x(:,1)+h*s(:,i),x(:,2)],dt,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,FCounter,VariableCSG) - Fx )/h;
    FCounter = FCounter + 1;
end


%Organize components for J
for i = 0:n-2
    J{1}(i+1) = s(i+1,mod(i+1,3)+1);
    J{2}(i+1) = s(i+1,mod(i,3)+1);
    J{3}(i+1) = s(i+2,mod(i,3)+1);
end
J{2}(n) = s(n,mod(n-1,3)+1);

end



function [h,Fh,FCounter] = LineSearch(h,dt,Fh,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,delh,setting,FCounter,VariableCSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function for lineseaching in the Newtons method.
%
% INPUTS:
% h           = the vector of all the pressure heads at each node at each 
%               time step.
% dt          = is the time vector,
%                   dt(1) is forward time step,
%                   dt(2) is backward time step.
% z           = the vector that contains all of the height of all the nodes.
% zAt51       = vector that contains infomation above the node closest to the
%               height z = 51,
%                   zAt51(1) is actual height at the node,
%                   zAt51(2) is the index of where the node is.
% Del         = the cell that contains all the differences,
%                   Del{1} is volume of each cell,
%                   Del{2} is distance between each node,
%                   del{3} is the time difference.
% R           = the rainfall on the suface (m)
% Kzz         = the vector that conatians the Kzz values for the different
%               sections in the model,
%                   Kzz(1) is aqufier,
%                   Kzz(2) is coal measure.
% Kxx         = the vector that conatians the Kxx values for the different
%               sections in the model,
%                   Kxx(1) is aqufier,
%                   Kxx(2) is coal measure.
% CSGon       = the state of CSG extraction, 1 is on and 0 is off.
% droughtOn   = the state of drought in the model, 1 is in drough and 0 is
%               out of drought.
% CSGDistance = the distance of CSG Plant.
% theta       = setting for the interation,
%                   theta(1) is for flux,
%                   theta(2) is for source.
% sigma       = the control on what aproximation will be used on the velocity
%               terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%               upwinding
% delh        = the full Newton step.
% setting     = string for what type of linesearching, takes:
%                   'simple' and 'two-point'
% FCounter    = counts of how many time F is computed.
% VariableCSG = the string that is the condition of whether the CSG plant
%               is varrying, takes:
%                   'on' and 'off'
% OUTPUTS:
% h           = the vector of all the pressure heads at each node at each 
%               time step.
% Fh          = F evaluted at h.
% FCounter    = counts of how many time F is computed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;                                      % initialise iterate
maxIt = 5;                                  % maximum iterations for the line search
lambda = 1;                                 % initialise lambda
hdag = h(:,1) + lambda*delh;                % the first h dagger 
[Fhdag,FCounter] = Ffunc([hdag,h(:,2)],dt,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,FCounter,VariableCSG); % F(h dagger)
sig0 = .1;                                  % control parameters for two point linesearching 
sig1 = .5;                                  % same ^
alpha = 10e-4;

switch setting
    case 'simple'
        while norm(Fh,2) <= norm(Fhdag,2) && i < maxIt
        i       = i+1;
        lambda  = 0.5*lambda; %Calculate new lambda
        
        hdag    = h(:,1) + lambda * delh; %Guess at new h
        [Fhdag,FCounter]   = Ffunc([hdag,h(:,2)],dt,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,FCounter,VariableCSG);
        end
        
    case 'two-point'
        while norm(Fh,2)*(1-alpha*lambda) <= norm(Fhdag,2) && i < maxIt
        i           = i+1;
        lambdaDagg  = (norm(Fh,2)*lambda)^2/(norm(Fhdag,2)^2 + ( 2*lambda - 1)*norm(Fh,2)^2); %Calculate new lambda
        
        if lambdaDagg < sig0*lambda %condistions on lambda so not to small or large
            lambda  = sig0 * lambda;
        elseif lambdaDagg > sig1*lambda
            lambda  = sig1 * lambda;
        else
            lambda  = lambdaDagg;
        end
        
        
        hdag = h(:,1) + lambda * delh; %Guess at new h
        [Fhdag,FCounter] = Ffunc([hdag,h(:,2)],dt,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,FCounter,VariableCSG);
        end
    otherwise
        error('Invalid String: Only "two-point" or "simple"')        
end

if i < maxIt 
    h = [hdag,h(:,2)];
    Fh = Fhdag; % parce the new values back if fell under iteration cap
else
    error('LineSearch:NoConverg','No Convergence under current timestep')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions for the Function Vector
function epsi = EpsiCal(h,dt,z,zAt51,Del,Kxx,CSGon,droughtOn,R,CSGDistance,VariableCSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the error in the water buget in the current time
% step.
%
% INPUTS:
% h           = the vector of all the pressure heads at each node at each 
%               time step.
% dt          = is the time vector,
%                   dt(1) is forward time step,
%                   dt(2) is backward time step.
% z           = the vector that contains all of the height of all the nodes.
% zAt51       = vector that contains infomation above the node closest to the
%               height z = 51,
%                   zAt51(1) is actual height at the node,
%                   zAt51(2) is the index of where the node is.
% Del         = the cell that contains all the differences,
%                   Del{1} is volume of each cell,
%                   Del{2} is distance between each node,
%                   del{3} is the time difference.
% Kxx         = the vector that conatians the Kxx values for the different
%               sections in the model,
%                   Kxx(1) is aqufier,
%                   Kxx(2) is coal measure.
% CSGon       = the state of CSG extraction, 1 is on and 0 is off.
% droughtOn   = the state of drought in the model, 1 is in drough and 0 is
%               out of drought.
% R           = the rainfall on the suface (m)
% CSGDistance = the distance of CSG Plant.
% VariableCSG = the string that is the condition of whether the CSG plant
%               is varrying, takes:
%                   'on' and 'off'
% OUTPUTS:
% epsi        = error in the water bugdet.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    epsi = NorthBoundaryFunc(R(1),h(size(h,1),1));
    for i = 1:length(z)
        epsi = epsi + Del{1}(i)*(phiFunc(h(i,1),z(i)) - phiFunc(h(i,2),z(i)))/Del{3} - (QFunc(1,dt,z(i),zAt51,h,Kxx,CSGon(1),droughtOn(1),R(1),phiFunc(h(i,1),z(i)),CSGDistance,VariableCSG));
    end
    
end

function [F,counter] = Ffunc(h,dt,z,zAt51,Del,R,Kzz,Kxx,CSGon,droughtOn,CSGDistance,theta,sigma,counter,VariableCSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function that is being solved in Newtons method.
%
% INPUTS:
% h           = the vector of all the pressure heads at each node at each 
%               time step.
% dt          = is the time vector,
%                   dt(1) is forward time step,
%                   dt(2) is backward time step.
% z           = the vector that contains all of the height of all the nodes.
% zAt51       = vector that contains infomation above the node closest to the
%               height z = 51,
%                   zAt51(1) is actual height at the node,
%                   zAt51(2) is the index of where the node is.
% Del         = the cell that contains all the differences,
%                   Del{1} is volume of each cell,
%                   Del{2} is distance between each node,
%                   del{3} is the time difference.
% R           = the rainfall on the suface (m)
% Kzz         = the vector that conatians the Kzz values for the different
%               sections in the model,
%                   Kzz(1) is aqufier,
%                   Kzz(2) is coal measure.
% Kxx         = the vector that conatians the Kxx values for the different
%               sections in the model,
%                   Kxx(1) is aqufier,
%                   Kxx(2) is coal measure.
% CSGon       = the state of CSG extraction, 1 is on and 0 is off.
% droughtOn   = the state of drought in the model, 1 is in drough and 0 is
%               out of drought.
% CSGDistance = the distance of CSG Plant.
% theta       = setting for the interation,
%                   theta(1) is for flux,
%                   theta(2) is for source.
% sigma       = the control on what aproximation will be used on the velocity
%               terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%               upwinding
% counter     = counts of how many time F is computed.
% VariableCSG = the string that is the condition of whether the CSG plant
%               is varrying, takes:
%                   'on' and 'off'
% OUTPUTS:
% F           = F evaluted at h.
% counter     = counts of how many time F is computed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Evaluates the vector function for the 1d dicretization at particular point

L       = length(h);
F       = zeros(L,1);



% Calculate the first valuse of the F(x) function
F(1) = phiFunc(h(1,1),z(1)) - phiFunc(h(1,2),z(1))...
        + Del{3} / Del{1}(1) * (theta(1) * ( qfunc(h(2,1),h(1,1),z(1),Del{2}(1),Kzz,sigma) - SouthBoundryFunc) ...
            + (1-theta(1)) *( qfunc(h(2,2),h(1,2),z(1),Del{2}(1),Kzz,sigma) - SouthBoundryFunc )  )...
        - Del{3} * (theta(2) * QFunc(1,dt,z(1),zAt51,h,Kxx,CSGon(1),droughtOn(1),R(1),phiFunc(h(1,1),z(1)),CSGDistance,VariableCSG) ...
            + (1 - theta(2)) * QFunc(2,dt,z(1),zAt51,h,Kxx,CSGon(2),droughtOn(2),R(2),phiFunc(h(1,2),z(1)),CSGDistance,VariableCSG));
    
% Calculate all but the last and the first valued of F(x)
for i = 2:L-1
    F(i) = phiFunc(h(i,1),z(i))- phiFunc(h(i,2),z(i)) ... 
        + Del{3} / Del{1}(i) * (  theta(1) * (qfunc(h(i+1,1),h(i,1),z(i),Del{2}(i),Kzz,sigma) - qfunc(h(i-1,1),h(i,1),z(i),-Del{2}(i-1),Kzz,sigma)) ...
            + (1-theta(1)) * (qfunc(h(i+1,2),h(i,2),z(i),Del{2}(i),Kzz,sigma) - qfunc(h(i-1,2),h(i,2),z(i),-Del{2}(i-1),Kzz,sigma))  ) ...
        - Del{3} * (theta(2) * QFunc(1,dt,z(i),zAt51,h,Kxx,CSGon(1),droughtOn(1),R(1),phiFunc(h(i,1),z(i)),CSGDistance,VariableCSG) ...
            + (1 - theta(2)) * QFunc(2,dt,z(i),zAt51,h,Kxx,CSGon(2),droughtOn(2),R(2),phiFunc(h(i,2),z(i)),CSGDistance,VariableCSG));    
end

% Calculate the last value of F(x)
F(L) =  phiFunc(h(L,1),z(L)) - phiFunc(h(L,2),z(L))...
        + Del{3} / Del{1}(L) * (  theta(1) * (NorthBoundaryFunc(R(1),h(L,1)) - qfunc(h(L-1,1),h(L,1),z(L),-Del{2}(L-1),Kzz,sigma)) ...
            + (1-theta(1)) * (NorthBoundaryFunc(R(2),h(L,2)) - qfunc(h(L-1,2),h(L,2),z(L),-Del{2}(L-1),Kzz,sigma))  ) ...
        - Del{3} * (theta(2) * QFunc(1,dt,z(L),zAt51,h,Kxx,CSGon(1),droughtOn(1),R(1),phiFunc(h(L,1),z(L)),CSGDistance,VariableCSG) ...
            + (1 - theta(2)) * QFunc(2,dt,z(L),zAt51,h,Kxx,CSGon(2),droughtOn(2),R(2),phiFunc(h(L,2),z(L)),CSGDistance,VariableCSG));
    
    counter = counter + 1;
end


function q = qfunc(h,hp,z,deltz,Kzz,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determines the flux between two nodes
%
% INPUTS:
% h     = the presure at next node.
% hp    = the presure at the focus node.
% z     = the height of the focus node.
% deltz = distance to next node.
% Kzz   = the vector that conatians the Kzz values for the different
%         sections in the model,
%           Kzz(1) is aqufier,
%           Kzz(2) is coal measure.
% sigma = the control on what aproximation will be used on the velocity
%          terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%          upwinding
% OUTPUTS:
% q = the flux between th nodes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if sigma == 2
    q = - harmmean([kFunc(h,z+deltz) kFunc(hp,z+deltz)]) * KzzFunc(Kzz,z,deltz,sigma,'harm') * ( (h-hp)/deltz + 1);
else
    if h < hp    
        q = - ( kFunc(h,z+deltz) - sigma*(kFunc(h,z+deltz) - kFunc(hp,z+deltz))/2) * KzzFunc(Kzz,z,deltz,sigma,'out') * ( (h-hp)/deltz + 1);
    else
        q = - ( kFunc(hp,z+deltz) - sigma*(kFunc(hp,z+deltz) - kFunc(h,z+deltz))/2) * KzzFunc(Kzz,z,deltz,sigma,'in') * ( (h-hp)/deltz + 1);
    end
end

end

function Kzz = KzzFunc(KzzVec,z,deltz,sigma,flow)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function determines what Kzz value is used in the flux calculation
%
% INPUTS:
% KzzVec = the vector that conatians the Kzz values for the different
%          sections in the model,
%            KzzVec(1) is aqufier,
%            KzzVec(2) is coal measure. 
% z      = the height of the node.
% deltz  = distance to next node.
% sigma  = the control on what aproximation will be used on the velocity
%          terms. 2 is harmonic averaging, 1 is averaging and, 0 is 
%          upwinding
% flow   = which direction the water is flowing, takes:
%           'in', 'out' and 'harm'
% OUTPUTS:
% Kzz    = Kzz to be used.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if deltz == 0
    error('Input Variable deltz, has been passed zero')
end


if z == 50
    switch flow
        case 'in'
            if deltz < 0
                Kzz = mean(KzzVec) - sigma * (mean(KzzVec) - KzzVec(2))/2;
            else
                Kzz = mean(KzzVec) - sigma * (mean(KzzVec) - KzzVec(1))/2;
            end
        case 'out'
            if deltz < 0
                Kzz = KzzVec(2) - sigma * (KzzVec(2) - mean(KzzVec))/2;
            else
                Kzz = KzzVec(1) - sigma * (KzzVec(1) - mean(KzzVec))/2;
            end
        case 'harm'
            if deltz < 0
                Kzz = harmmean([KzzVec(2) mean(KzzVec)]);
            else
                Kzz = harmmean([KzzVec(1) mean(KzzVec)]);
            end
        otherwise
            error('Input Variable flow, has been passed some thing other than "in" or "out"')
    end
elseif z + deltz < 50    
    Kzz = KzzVec(2);
else
    Kzz = KzzVec(1);    
end

end

function Q = QFunc(tStep,dt,z,zAt51,h,Kxx,CSGon,droughtOn,R,phi,DeltXc,VariableCSG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is the fuction for the source term
%
% INPUTS:
% tStep       = timestep, 1 is forward and 0 is backward.
% dt          = is the time vector,
%                 dt(1) is forward time step,
%                 dt(2) is backward time step.
% z           = the height of the node.
% zAt51       = vector that contains infomation above the node closest to the
%               height z = 51,
%                 zAt51(1) is actual height at the node,
%                 zAt51(2) is the index of where the node is.
% h           = the vector of all the pressure heads at each node at each 
%               time step.
% Kxx         = the vector that conatians the Kxx values for the different
%               sections in the model,
%                 Kxx(1) is aqufier,
%                 Kxx(2) is coal measure.
% CSGon       = the state of CSG extraction, 1 is on and 0 is off.
% droughtOn   = the state of drought in the model, 1 is in drough and 0 is
%               out of drought.
% R           = the rainfall on the suface (m)
% phi         = the water content of the soil at the node.
% DeltXc      = the distance of CSG Plant
% VariableCSG = the string that is the condition of whether the CSG plant
%               is varrying, takes:
%                   'on' and 'off'
% OUTPUTS:
% Q           = the source at the node
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ha          = (71.5-70.5)/2*cos(dt(tStep)/6*pi)+(71.5+70.5)/2;  % head of the river 
DeltXa      = 1000;                                             % distance for river
l           = 10;                                               % EVT region
L           = 80;                                               % height of the z region                                                         

% CSG head
switch VariableCSG 
    case 'on'                                                                   
        if dt(tStep) < 5*365           
            Hc          = 20;
            %disp('The CSG is running with head pressure 20.')
        else
            Hc          = 10; 
            %disp('The CSG is running with head pressure 10.')
        end
    case 'off'
        Hc = 20; 
end


% This is for CSG extraction.
if z == zAt51(2)
    % Source for the aquifer.
    Q = (1-droughtOn)*Kxx(1) * (Ha - h(zAt51(1),tStep) - zAt51(2)) / DeltXa ;
elseif z == 0
    % Source for the CSG plant.
    Q = CSGon * Kxx(2) * (Hc - h(1,tStep) - 0) / DeltXc ;
else 
    Q = 0;
end

% This is for the EVT effects.
if z >= L-l && phi/0.33 > 0.5
    Q = Q+0.1*R*(z-L+l)^2/l^2;
end
   

end

%% Boundary Condition Function

function qs = SouthBoundryFunc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function that controls the flux at the south boundary
%
% INPUTS:
% N\A
% OUTPUTS:
% qs = the flux at the south boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    qs = 0;
end

function qn = NorthBoundaryFunc(R,h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the function that controls the flux at the north boundary
%
% INPUTS:
% R  = the rainfall on the suface (m)
% h  = the pressure at the north boundary node
% OUTPUTS:
% qn = the flux at the north boundary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Shock = @(x) 1/(1+exp(35*x));
    qn = Shock(h) * R;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Base Functions

function phi = phiFunc(h,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The phi function for the model.
%
% INPUTS:
% h = pressure at the node.
% z = height of the node.
% OUTPUTS:
% phi = phi evaluated at the node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Moisture content
sat = ConstChooser('sat',z);
res = ConstChooser('res',z);

if h < 0
    phi = res + SFunc(h,z) * (sat - res);
else
    phi = sat;
end

end

function k = kFunc(h,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The k function for the model.
%
% INPUTS:
% h = pressure at the node.
% z = height of the node.
% OUTPUTS:
% k = k evaluated at the node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = ConstChooser('n',z);
m = 1 - 1/n;

if h < 0
    k = sqrt(SFunc(h,z)) * (1 - (1 - SFunc(h,z)^(1/m) )^m )^2;
else
    k = 1;
end

end

function S = SFunc(h,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The S function for the model.
%
% INPUTS:
% h = pressure at the node.
% z = height of the node.
% OUTPUTS:
% S = S evaluated at the node.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha = ConstChooser('alpha',z);
n = ConstChooser('n',z);
m = 1 - 1/n;

if h < 0
    S  = (1 + (-alpha*h)^n)^-m;
else
    S = 1;
end

end

function x = ConstChooser(string,z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose the correct constant for a node in the model
%
% INPUTS:
% z      = the height of the node
% string = what constant is needed, takes:
%           'sat', 'res', 'alpha', and 'n'
% OUTPUTS:
% x      = the constant need.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch string
    case 'sat'
        if z < 50
            x = .1;
        else
            x = .33;
        end
    case 'res'
        if z < 50
            x = 0;
        else
            x = 0;
        end
    case 'alpha'
        if z < 50
            x = 1.43;
        else
            x = 1.43;
        end
    case 'n'
        if z < 50
            x = 1.51;
        else
            x = 1.51;
        end
    otherwise
        error('Incorrect String: Only can input sat, res, alpha, n.')
end
end

function sourceStates = DeterminState(timeState)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function deterimes whether the function is going into, going out,
% fully in or fully out and set the source terms accordingly
%
% INPUTS:
% timeState    = the string that tells what state the models in, can take: 
%                   'FullyIn', 'GoingOut', 'FullyOut', and 'GoingIn'
% OUTPUTS:
% sourceStates = a vector that contains the states that the sources terms 
%                are in,
%                   sourceStates(1) is forward time step,
%                   sourceStates(2) is backward time step.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch timeState
    case 'FullyIn'
        sourceStates = [1;1];
    case 'GoingOut'
        sourceStates = [0;1];
    case 'FullyOut' 
        sourceStates = [0;0];
    case 'GoingIn'
        sourceStates = [1;0];
    otherwise
        eString = strcat('String "',timeState,'" is not accepted');
        error(eString)
end
end
