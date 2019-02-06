%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MXB324 
% Group 1: Johnathan Adams, James Beattie, Georgia Harman, Thyraphol Sutanujinda
% Dalby Groundwater model 
% Main Script Operations Script
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the main operations script for Group 1's groundwater model solver. 
% From this script the user can change the parameters of the model and also
% request for various parts in the modelling process to be visualised or
% not. The operations script can be broken up into five parts,
% (1) Parameter Initialisation
% (2) Tensor Value Estimation and 2D Grid Generation
% (3) 1D Grid Generation
% (4) Pressure Head Solver
% (5) Statistics and Visualisation of Solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (1) Parameter Initialisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
% Grid Parameters

% One Dimensionsal Grid
n1d     = 160;                              % resolution of the z grid.
L       = 80;                               % length of the one d grid.
type_1d = 'uniform';                        % the type of grid to be generated 
viz_1d  = 'off';                            % 1D visualisation

% Two Dimensional Grid
m2d             = 21;                       % number of nodes in a unit cell in the x direction 
n2d             = 21;                       % number of nodes in a unit cell in the z direction
tol             = 10^-8;                    % the tolerance for the 2D GMRES
maxiters        = 300;                      % the maximum number of iterations for the 2D GMRES
precondition    = 'jacobi';                 % the preconditioning matrix
viz_2d          = 'off';                    % 2D visualisation
type_2d         = 'uniform';                % the type of grid to be generated for the 2d problem
parms_2d        = {m2d,n2d,tol,maxiters,precondition,viz_2d,type_2d};


% Source Term Paramaters
CSGPeriods     = [0 0];
droughtPeriods = [0 0];
CSGDistance    = 5000;
rainSetting    = 'High Frequency';
VariableCSG    = 'on';

% Time Parameters
tInterval   = [0; 50*365];  %[Starting time (days) ; Finishing time (days)]
delt        = 0.365;        %Initial Time Step (days)

% Newton's Method Parameters     
solveSetting    = 'TDMASolve';
m               = 3;              % Full Newton method = 1, not equal to one, Shamanski
theta           = [1;1];            % theta=1 backward euler, theta=1/2 CN
sigma           = 1;                % sigma = 2 is harmonic averaging, sigma = 1 is averaging, sigma = 0 is upwinding

% Analysis Parameters
type    = 'Statistics';
choice  = 'PressureHead';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (2) Tensor Value Estimates: The 2D Problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Kxx,Kzz]  = ExtractTensorValues(parms_2d);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (3) 1D Grid Generator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Delta,z,hInitial]   = OneDGridGen(delt,n1d,L,viz_1d,type_1d);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (4) Newton Solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[hStorage,runTime,timeSteps,FEvals,espis] = ...
    OneDModel(tInterval,z,Delta,Kzz,Kxx,hInitial,CSGPeriods,droughtPeriods,CSGDistance,rainSetting,solveSetting,VariableCSG,m,sigma,theta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% (5) Point Statistics and Visualisation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

StatisticsandVisualisation(hStorage,runTime,timeSteps,FEvals,tInterval,z,espis,type,choice)
StatisticsandVisualisation(hStorage,runTime,timeSteps,FEvals,tInterval,z,espis,'Visualisation','PressureHead')

