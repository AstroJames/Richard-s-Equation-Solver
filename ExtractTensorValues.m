function [Kxx,Kzz]= ExtractTensorValues(parms)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExtractTensorValues calculates the effetive conductivities for both
% alluvium and walloon coal. The function generates a 2D grid, constructs
% linear systems, stores the linear systems, solves the linear systems, and 
% performs a bilinear interpolation. This performed by calling 5 other 
% functions to perform this:twod_gridgenerator, construct_A_b, 
% twod_storage,twod_solver and tensor_solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% parms         =   a cell containing an input variable to the 2D problem,
%                   where:
%                   -   parms{1} contain a scalar m2d, which is the number
%                   of nodes in a unit cell in the x direction,
%                   -   parms{2} contain a scalar n2d, which is the number
%                   of nodes in a unit cell in the z direction,
%                   -   parms{3} contains a scalar tol, which is the
%                   tolerance for the 2D GMRES,
%                   -   parms{4} contain a scalar maxiters, which is the
%                   maximum number of iterations for the 2D GMRES,
%                   -   parms{5} contain a string precondition, which is
%                   the type of preconditioning used in 2D GMRES,
%                   -   parms{6} contain a string viz_2d, which is the
%                   status of 2D visualisation
%                   -   parms{7} contain a string type_2d, which is the
%                   type of grid to be generated for 2D problem.
% Output:
% Kxx           =   a vector [Kxx_a,Kxx_c] where,
%                   -   Kxx_a is the effective conductivity of alluvium
%                   aquifer in the x-direction,
%                   -   Kxx_c is the effective conductivity of walloon coal
%                   in the x-direction.
% Kzz           =   a vector [Kxx_a,Kxx_c] where,
%                   -   Kzz_a is the effective conductivity of alluvium
%                   aquifer in the z-direction,
%                   -   Kzz_c is the effective conductivity of walloon coal
%                   in the z-direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract input
m2d             = parms{1};
n2d             = parms{2};
tol             = parms{3};
maxiters        = parms{4};
precondition    = parms{5};
viz_2d          = parms{6};
type_2d         = parms{7};

disp('BEGINNING TO EXTRACT TENSOR VALUES')
fprintf('You have selected the preconditioning method %s.\n',precondition)


%% 2D Grid Generator
fprintf('Generating %f by 1%f grid\n',m2d,n2d)
[x2d,z2d,del,m2d,n2d]       = twod_gridgenerator(m2d, n2d, type_2d, viz_2d);
x0                          = zeros((m2d-1)*(n2d-1),1);  % the initial guess for the 2d GMRES

% Alluvium Tensor Values
disp('Solving for the alluvium tensor values')
tic
[A,b1,b2,r]                 = Construct_A_b(x2d,z2d,del,'a'); % Construct linear system
[A, b1, b2, p_inverse]      = twod_storage(A,b1, b2); % Convert storage type
[u1_a,u2_a,residual,k]      = twod_solver(x0,A,b1,b2,tol,precondition,maxiters,m2d,n2d,p_inverse,r); % Solve for u1 and u2
[Kxx_a,Kzz_a]               = Tensor_solver(u1_a,u2_a,x2d,z2d,del,'a'); % Calculate effective conductivity terms
toc
fprintf('The alluvium aquifer effective conductivity are %f in the x direction and %f in the z direction\n',Kxx_a,Kzz_a)
disp('Complete.')
% Walloon Coal Tensor Values
disp('Solving for the coal measure tensor values')
tic
[A,b1,b2,r]                   = Construct_A_b(x2d,z2d,del,'c'); % Construct linear system
[A, b1, b2, p_inverse]      = twod_storage(A,b1, b2); % Convert storage type
[u1_c,u2_c,residual,k]      = twod_solver(x0,A,b1,b2,tol,precondition,maxiters,m2d,n2d, p_inverse,r); % Solve for u1 and u2
[Kxx_c,Kzz_c]               = Tensor_solver(u1_c,u2_c,x2d,z2d,del,'c'); % Calculate effective conductivity terms
toc
fprintf('The walloon coal measure effective conductivity are %f in the x direction and %f in the z direction\n',Kxx_c,Kzz_c)
disp('Complete.')

% Visualisation of u1 and u2
switch viz_2d
    case 'on'
        disp('You have select visualisation on for u1 and u2 visualisation.')
        
        % Alluvium Aquifier
        u_visualisation(x2d,z2d,u1_a,u2_a)
        subplot(1,2,1)
        title('Plot of the solution for u1 across the Alluvium Aquifier Unit Cell')
        subplot(1,2,2)
        title('Plot of the solution for u2 across the Alluvium Aquifier Unit Cell')
        
        % Walloon Coal
        u_visualisation(x2d,z2d,u1_c,u2_c)
        subplot(1,2,1)
        title('Plot of the solution for u1 across the Walloon Coal Unit Cell')
        subplot(1,2,2)
        title('Plot of the solution for u2 across the Walloon Coal Unit Cell')
        
    case 'off'
        disp('You have select visualisation off for u1 and u2 visualisation.')
end

% Store all Tensor Values
Kxx = [Kxx_a,Kxx_c];
Kzz = [Kzz_a,Kzz_c];
disp('TENSOR VALUE EXTRACT COMPLETE')


end