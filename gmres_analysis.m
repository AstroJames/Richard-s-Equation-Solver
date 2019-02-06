% GMRES ANALYSIS
% this script performs the GMRES precondition analysis outlined in the final report
% the analysis determines, run times, number of iterations and residuals
% additionally the optimal omega and 


% Two Dimensional Grid Parameters
m2d             = 21;               % number of nodes in a unit cell in the x direction 
n2d             = 21;               % number of nodes in a unit cell in the z direction
tol             = 10^-8;            % the tolerance for the 
maxiters        = 200;              % the maximum number of iterations for the 2D GMRES 
precondition    = 'incomplete';         % the preconditioning matrix
viz_2d          = 'off';            % 2D visualisation
type_2d         = 'uniform';        % the type of grid to be generated for the 2d problem
%%
%%%%%%%%%% Run Times, Iterations and Residuals %%%%%%%%%%%%%%%%%%%%%%


% set up grid
[x2d,z2d,del,m2d,n2d]       = twod_gridgenerator(m2d, n2d, type_2d, viz_2d);
x0                          = zeros((m2d-1)*(n2d-1),1);  % the initial guess for the 2d GMRES

%%%%%%%%%%%%%%%% ALLUVIUM %%%%%%%%%%%%%%%%%%%%%%%
% set up linear system
[A1,b1_1,b2_1,r]                  = Construct_A_b(x2d,z2d,del,'a');
[A1, b1_1, b2_1, p_inverse]     = twod_storage(A1,b1_1, b2_1);
tic; % record GMRES time
[u1_a,u2_a,residual,k]          = twod_solver(x0,A1,b1_1,b2_1,tol,precondition,maxiters,m2d,n2d, p_inverse,r);
precondition_toc_a = toc;
  

% display statistics
disp(['The GMRES run time for solving the Alluvium linear system using the ' precondition ' precondition is ' num2str(precondition_toc_a)])
disp(['The number of iterations solving the Alluvium linear system using the ' precondition ' precondition is ' num2str(k(1))])
disp(['The GMRES residual for solving the Alluvium linear system using the ' precondition ' precondition is ' num2str(residual(1))])


%%%%%%%%%%%%%%%%%% WALLOON COAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up linear system
[A2,b1_2,b2_2,r]                   = Construct_A_b(x2d,z2d,del,'c');
[A2, b1_2, b2_2, p_inverse] = twod_storage(A2,b1_2, b2_2);
tic; %record GMRES time
[u1_c,u2_c,residual,k_b] = twod_solver(x0,A2,b1_2,b2_2,tol,precondition,maxiters,m2d,n2d,p_inverse, r);
precondition_toc_b = toc;

%display statistics
disp(['The GMRES run time for solving the Walloon linear system using the ' precondition ' precondition is ' num2str(precondition_toc_a)])
disp(['The number of iterations solving the Walloon linear system using the ' precondition ' precondition is ' num2str(k(2))])
disp(['The GMRES residual for solving the Walloon linear system using the ' precondition ' precondition is ' num2str(residual(2))])




%%
%%%%%%%%%%%%%%%%%%%%% CALCULATE OPTIMAL OMEGA FOR SOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% alluvium aquifer
[A1,b1_1,b2_1]                  = Construct_A_b(x2d,z2d,del,'a'); % construct linear system
J_a                             = diag(diag(A1)); % set up Jacobian precondition
omega_a                         = 2/ (1 + sqrt(max(abs(eig(J_a)))^2)) ; % omega caluclation



[A2,b1_2,b2_2]                   = Construct_A_b(x2d,z2d,del,'c'); %construct linear system
J_c                              = diag(diag(A2)); %set up jacobian precondition
omega_c                          = 2/ (1 + sqrt(max(abs(eig(J_c)))^2));

% display omega
disp(['The optimal omega for Alluvium is ' num2str(omega_a)])
disp(['The optimal omega for Walloon is ' num2str(omega_c)])

