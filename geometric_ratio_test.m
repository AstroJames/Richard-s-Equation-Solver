%% Geometric spacing ratio's convergence test (Analysis purpose only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script file perform an analysis on effect of geometric spacing ratio
% on the convergence of the 2D problem. Initially, it solves for the
% uniform grid problem with 41x41 and 61x61 grids as a controlled test.
% Then, it performs the same test on the bilateral geometric grid of the
% same size, varying the spacing ratio in each iteration. The relevant
% result are then plot against the spacing ratio and compared to uniform
% grid problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: since special bilateral geometric grid is designed specifically
% for Group 1 alluvium aquifer, only tests which involves alluvium aquifer
% will have 6 plots instead of 4 on the same figure, and 2 of them may not
% be valid for unit cell geometry.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WARNING: Code runtime will take several hours to run.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Clear MATLAB
close all
clear all
clc

%% Set up parameters
max_iter        = 50;
iter_length     = 1:max_iter;
type            = 'bilateral';
viz             = 'off';
Kxx_a_length    = zeros(max_iter,2);
Kxx_c_length    = zeros(max_iter,2);
Kzz_a_length    = zeros(max_iter,2);
Kzz_c_length    = zeros(max_iter,2);
Kxx_as_length   = zeros(max_iter,2);
Kzz_as_length   = zeros(max_iter,2);

%% Uniform comparison
% 41 x 41 grid size
[x,z,del,m,~]       = twod_gridgenerator(41, 41,'uniform','off');
% Alluvium Aquifer
[A,b1,b2]           = Construct_A_b(x,z,del,'a');
[u1,u2]             = twod_analysis(A,b1,b2,m,m);
[Kxx_au,Kzz_au]     = Tensor_solver(u1,u2,x,z,del,'a');
% Walloon Coal
[A,b1,b2]           = Construct_A_b(x,z,del,'c');
[u1,u2]             = twod_analysis(A,b1,b2,m,m);
[Kxx_cu,Kzz_cu]     = Tensor_solver(u1,u2,x,z,del,'c');

% 61 x 61 grid size
[x,z,del,m,~]       = twod_gridgenerator(61, 61,'uniform','off');
% Alluviun Aquifer
[A,b1,b2]           = Construct_A_b(x,z,del,'a');
[u1,u2]             = twod_analysis(A,b1,b2,m,m);
[Kxx_au2,Kzz_au2]   = Tensor_solver(u1,u2,x,z,del,'a');
% Walloon Coal
[A,b1,b2]           = Construct_A_b(x,z,del,'c');
[u1,u2]             = twod_analysis(A,b1,b2,m,m);
[Kxx_cu2,Kzz_cu2]   = Tensor_solver(u1,u2,x,z,del,'c');

%% Bilateral geometric grid test
for i = 1:max_iter
    % Set up spacing ratio
    s = 1+i/max_iter;
    
    % Standard bilateral geometric grid
    % 41x41 Grid
    [x,z,del,~,m]   = twod_gridgenerator(41,41, type, viz, s);
    % Alluvium Aquifer
    [A,b1,b2]       = Construct_A_b(x,z,del,'a');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_a,Kzz_a]   = Tensor_solver(u1,u2,x,z,del,'a');
    % Walloon Coal
    [A,b1,b2]       = Construct_A_b(x,z,del,'c');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_c,Kzz_c]   = Tensor_solver(u1,u2,x,z,del,'c');
    
    % 61x61 Grid
    [x,z,del,~,m]   = twod_gridgenerator(61,61, type, viz, s);
    % Alluvium Aquifer
    [A,b1,b2]       = Construct_A_b(x,z,del,'a');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_ar,Kzz_ar] = Tensor_solver(u1,u2,x,z,del,'a');
    % Walloon Coal
    [A,b1,b2]       = Construct_A_b(x,z,del,'c');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_cr,Kzz_cr] = Tensor_solver(u1,u2,x,z,del,'c');
    
    % Special Bilateral Geometric grid
    % 41x41 grid
    [x,z,del,~,m]   = twod_gridgenerator(41,41, type, viz, 'a', s);
    % Alluvium Aquifer
    [A,b1,b2]       = Construct_A_b(x,z,del,'a');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_as,Kzz_as]   = Tensor_solver(u1,u2,x,z,del,'a');
    
    % 61x61 grid
    [x,z,del,~,m]   = twod_gridgenerator(61,61, type, viz, 'a', s);
    % Alluvium Aquifer
    [A,b1,b2]       = Construct_A_b(x,z,del,'a');
    [u1,u2]         = twod_analysis(A,b1,b2,m,m);
    [Kxx_asr,Kzz_asr]   = Tensor_solver(u1,u2,x,z,del,'a');
    
    % Store Data
    sprintf('iteration number %d', i)
    Kxx_a_length(i,:)   = [Kxx_a,Kxx_ar];
    Kxx_c_length(i,:)   = [Kxx_c,Kxx_cr];
    Kzz_a_length(i,:)   = [Kzz_a,Kzz_ar];
    Kzz_c_length(i,:)   = [Kzz_c,Kzz_cr];
    Kxx_as_length(i,:)  = [Kxx_as,Kxx_asr];
    Kzz_as_length(i,:)  = [Kzz_as,Kzz_asr];
end

%% Visualisation
% Set up x-axis
ratio_length = 1+iter_length/max_iter;

% Kxx_a convergence plots
figure
plot(ratio_length,Kxx_a_length(:,1))
hold on
plot(ratio_length,Kxx_a_length(:,2))
plot(ratio_length,Kxx_as_length(:,1))
plot(ratio_length,Kxx_as_length(:,2))
plot(ratio_length([1,end]),[Kxx_au,Kxx_au])
plot(ratio_length([1,end]),[Kxx_au2,Kxx_au2])
title('Convergence Test: Alluvium Aquifer Kxx')
xlabel('Geometric Ratio')
ylabel('Effective Kxx')
legend('41x41 geometric','61x61 geometric','41x41 special geometric','61x61 special geometric','41x41 uniform','61x61 uniform')

% Kzz_a convergence plots
figure
plot(ratio_length,Kzz_a_length(:,1))
hold on
plot(ratio_length,Kzz_a_length(:,2))
plot(ratio_length,Kzz_as_length(:,1))
plot(ratio_length,Kzz_as_length(:,2))
plot(ratio_length([1,end]),[Kzz_au,Kzz_au])
plot(ratio_length([1,end]),[Kzz_au2,Kzz_au2])
title('Convergence Test: Alluvium Aquifer Kzz')
xlabel('Geometric Ratio')
ylabel('Effective Kxx')
legend('41x41 geometric','61x61 geometric','41x41 special geometric','61x61 special geometric','41x41 uniform','61x61 uniform')

% Kxx_c convergence plots
figure
plot(ratio_length,Kxx_c_length(:,1))
hold on
plot(ratio_length,Kxx_c_length(:,2))
plot(ratio_length([1,end]),[Kxx_cu,Kxx_cu])
plot(ratio_length([1,end]),[Kxx_cu2,Kxx_cu2])
title('Convergence Test: Walloon Coal Kxx')
xlabel('Geometric Ratio')
ylabel('Effective Kxx')
legend('41x41 geometric','61x61 geometric','41x41 uniform','61x61 uniform')

% Kzz_c convergence plots
figure
plot(ratio_length,Kzz_c_length(:,1))
hold on
plot(ratio_length,Kzz_c_length(:,2))
plot(ratio_length([1,end]),[Kzz_cu,Kzz_cu])
plot(ratio_length([1,end]),[Kzz_cu2,Kzz_cu2])
title('Convergence Test: Walloon Coal Kzz')
xlabel('Geometric Ratio')
ylabel('Effective Kxx')
legend('41x41 geometric','61x61 geometric','41x41 uniform','61x61 uniform')
