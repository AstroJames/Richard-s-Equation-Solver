%% Grid types convergence test (Analysis purpose only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script file perform an analysis on the convergence rate of different
% grid type for the 2D problem. For the different grid type, this function
% solves the effective conductivity terms for alluvium aquifer and walloon
% coal then repeats the process on a finer grid. The result are the plotted
% according to the step-size used to set up the grid. For consistency
% purpose, the number of nodes in the x and z-direction are set to be
% equal, and the grid ratio used for bilateral geometric grid are the
% approximated optimal value for the relevant test.
% The result of the analysis are split into 4 parts:
% 1.) Check the convergence of the uniform grid problem. This result will
% is used as the controlled test.
% 2.) Compare the convergence of the random grid problem to the controlled
% test.
% 3.) Compare the convergence of the optimised bilateral geometric grid
% problem to the controlled test.
% 4.) Compare the convergence of the optimised special bilateral geometric
% grid problem for alluvium aquifer to the controlled test.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Code runtime may take a few hours to run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear MATLAB
close all
clear all
clc
%% Convergence tests

% Set up parameters
max_iter            = 15;
iter_length         = 1:max_iter;
Kxx_length          = zeros(max_iter,2); % uniform Kxx
Kzz_length          = zeros(max_iter,2); % uniform Kzz 
Kxx_length_r        = zeros(max_iter,2); % random Kxx
Kzz_length_r        = zeros(max_iter,2); % random Kzz
Kxx_length_b        = zeros(max_iter,2); % bilateral Kxx
Kzz_length_b        = zeros(max_iter,2); % bilateral Kzz
K_length_bs         = zeros(max_iter,2); % special bilateral
viz                 = 'off';

for i = 1:max_iter
    % Set up grid resolution %
    m               = 4*i+1;
    % Uniform Grid %
    type            = 'uniform';
    [x,z,del,~,m]   = twod_gridgenerator(m,m, type, viz);
    % Alluvium Aquifer
    [A,b1,b2]       = Construct_A_b(x,z,del,'a');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_a,Kzz_a]   = Tensor_solver(u1,u2,x,z,del,'a');
    % Walloon Coal
    [A,b1,b2]       = Construct_A_b(x,z,del,'c');
    [u1,u2]       	= twod_analysis(A,b1,b2,m,m);
    [Kxx_c,Kzz_c]  	= Tensor_solver(u1,u2,x,z,del,'c');
    
    % Random Grid %
    type            = 'random';
    [x,z,del,~,m]   = twod_gridgenerator(m,m, type, viz);
    % Alluvium Aquifer
    [A,b1,b2]       = Construct_A_b(x,z,del,'a');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_ar,Kzz_ar] = Tensor_solver(u1,u2,x,z,del,'a');
    % Walloon Coal
    [A,b1,b2]       = Construct_A_b(x,z,del,'c');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_cr,Kzz_cr] = Tensor_solver(u1,u2,x,z,del,'c');
    
    % Bilateral Geometric grid %
    type            = 'bilateral';
    % Alluvium Aquifer
    [x,z,del,~,m]       = twod_gridgenerator(m,m, type, viz, 1.2820);
    [A,b1,b2]           = Construct_A_b(x,z,del,'a');
    [u1,u2]             = twod_analysis(A,b1,b2,m,m);
    [Kxx_ab,Kzz_ab]     = Tensor_solver(u1,u2,x,z,del,'a');
    % Walloon Coal
    [x,z,del,~,m]       = twod_gridgenerator(m,m, type, viz, 1.6141);
    [A,b1,b2]           = Construct_A_b(x,z,del,'c');
    [u1,u2]             = twod_analysis(A,b1,b2,m,m);
    [Kxx_cb,Kzz_cb]     = Tensor_solver(u1,u2,x,z,del,'c');
    
    % Special Bilateral Geometric grid (Alluvium Aquifer only) %
    [x,z,del,~,m]       = twod_gridgenerator(m,m, type, viz,'a', 1.1960);
    [A,b1,b2]           = Construct_A_b(x,z,del,'a');
    [u1,u2]             = twod_analysis(A,b1,b2,m,m);
    [Kxx_abs,Kzz_abs]   = Tensor_solver(u1,u2,x,z,del,'a');
    
    % Store data %
    sprintf('iteration number %d', i)
    Kxx_length(i,:)     = [Kxx_a,Kxx_c];
    Kzz_length(i,:)     = [Kzz_a,Kzz_c];
    Kxx_length_r(i,:)   = [Kxx_ar,Kxx_cr];
    Kzz_length_r(i,:)   = [Kzz_ar,Kzz_cr];
    Kxx_length_b(i,:)   = [Kxx_ab,Kxx_cb];
    Kzz_length_b(i,:)   = [Kzz_ab,Kzz_cb];
    K_length_bs(i,:)    = [Kxx_abs,Kzz_abs];
end

%% Uniform convergence plot
figure
subplot(2,2,1)
plot(iter_length,Kxx_length(:,1))
title('Convergence Test: Alluvium Aquifer (Kxx)')
xlabel('step size')
ylabel('Effective K')

subplot(2,2,2)
plot(iter_length,Kzz_length(:,1))
title('Convergence Test: Alluvium Aquifer (Kzz)')
xlabel('step size')
ylabel('Effective K')

subplot(2,2,3)
plot(iter_length,Kxx_length(:,2))
title('Convergence Test: Walloon Coal (Kxx)')
xlabel('step size')
ylabel('Effective K')

subplot(2,2,4)
plot(iter_length,Kzz_length(:,2))
title('Convergence Test: Walloon Coal (Kzz)')
xlabel('step size')
ylabel('Effective K')

%% Random convergence plot

figure
subplot(2,2,1)
plot(iter_length,Kxx_length(:,1))
hold on
plot(iter_length,Kxx_length_r(:,1))
title('Convergence Test: Alluvium Aquifer (Kxx)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Random')

subplot(2,2,2)
plot(iter_length,Kzz_length(:,1))
hold on
plot(iter_length,Kzz_length_r(:,1))
title('Convergence Test: Alluvium Aquifer (Kzz)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Random')

subplot(2,2,3)
plot(iter_length,Kxx_length(:,2))
hold on
plot(iter_length,Kxx_length_r(:,2))
title('Convergence Test: Walloon Coal (Kxx)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Random')

subplot(2,2,4)
plot(iter_length,Kzz_length(:,2))
hold on
plot(iter_length,Kzz_length_r(:,2))
title('Convergence Test: Walloon Coal (Kzz)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Random')

%% Bilateral Geometric convergence plot

figure
subplot(2,2,1)
plot(iter_length,Kxx_length(:,1))
hold on
plot(iter_length,Kxx_length_b(:,1))
title('Convergence Test: Alluvium Aquifer (Kxx)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Bilateral Geometric')

subplot(2,2,2)
plot(iter_length,Kzz_length(:,1))
hold on
plot(iter_length,Kzz_length_b(:,1))
title('Convergence Test: Alluvium Aquifer (Kzz)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Bilateral Geometric')

subplot(2,2,3)
plot(iter_length,Kxx_length(:,2))
hold on
plot(iter_length,Kxx_length_b(:,2))
title('Convergence Test: Walloon Coal (Kxx)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Bilateral Geometric')

subplot(2,2,4)
plot(iter_length,Kzz_length(:,2))
hold on
plot(iter_length,Kzz_length_b(:,2))
title('Convergence Test: Walloon Coal (Kzz)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform','Bilateral Geometric')

%% Special Bilateral Geometric convergence plot
figure
subplot(1,2,1)
plot(iter_length,Kxx_length(:,1))
hold on
plot(iter_length,Kxx_length_b(:,1))
plot(iter_length,K_length_bs(:,1))
title('Convergence Test: Alluvium Aquifer (Kxx)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform', 'Bilateral Geometric','Bilateral Geometric adjusted')

subplot(1,2,2)
plot(iter_length,Kzz_length(:,1))
hold on
plot(iter_length,Kzz_length_b(:,1))
plot(iter_length,K_length_bs(:,2))
title('Convergence Test: Alluvium Aquifer (Kzz)')
xlabel('step size')
ylabel('Effective K')
legend('Uniform', 'Bilateral Geometric','Bilateral Geometric adjusted')