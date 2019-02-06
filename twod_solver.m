function [u1,u2,residual,k] = twod_solver(x0,A,b1,b2,tol,precondition,maxiters,m,n,p_inverse,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% twod_solver solves the linear system Au = b for u1 and u2, then reordere
% and restore the solution u in a matrix form. Preconditioned GMRES is used
% to iteratively solve the linear system. Number of iterations and residual
% are record and output at the end of the function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% x0            =   initial guess vector of x the same size as b.
% A             =   a coefficnet matrix A stored using RCM ordering as a
%                   sparse matrix
% b1            =   the b vector for j = 1 with RCM ordering.
% b2            =   the b vector for j = 2 with RCM ordering.
% tol           =   tolerance level for 2D GMRES solver
% precondition  =   precondition matrix for 2D GMRES solver, where
%                   - 'no_pre' use no precondition,
%                   - 'jacobi' use Jacobi precondition,
%                   - 'gauss-siedel use Gauss-Sidel precondition,
%                   - 'SOR' use SOR precondition,
%                   - 'incomplete' use incomplete LU factorisation
%                   precondition.
% maxiters      =   the maximum number of iterations.
% m             =   the number of nodes in the x-direction.
% n             =   The number of nodes in the z-direction.
% p_inverse     =   the inverse RCM ordering permutation vector.
% r             =   the type of the material of interest, where
%                   - 'a' is the Alluvium Aquifer,
%                   - 'c' is the Walloon Coal Measure,
%                   - 't' is the benchmark test case.
% Output:
% u1        =   the numerical solution to the 2D PDE problem in
%               x-direction, stored in an m x n matrix.
% u2        =   the numerical solution to the 2D PDE problem in
%               z-direction, stored in an m x n matrix.
% residual  =   a vector containing two residuals of the iterative solution
%               for u1 and u2 respectively.
% k         =   a vector containing two final iteration counts for u1 and
%               u2 respectively.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[u1,residual1, k1]    = precondition_gmres(x0, A, b1, tol, precondition,maxiters,r);
[u2,residual2, k2]    = precondition_gmres(x0, A, b2, tol, precondition,maxiters,r);
u1 = u1(p_inverse);
u2 = u2(p_inverse);

u1 = vect_to_mat(u1,m,n);
u2 = vect_to_mat(u2,m,n);
residual = [residual1,residual2];
k = [k1,k2];

end

function [u_new] = vect_to_mat(u_old,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% vect_to_mat converts the storage of the solution to the homogenisation
% problem from its vector form to a matrix form. This allow for simpler
% indentation when solving for effective conductivities and during
% visualisation of the solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% u_old     =   the solution of the homogenisation problem stored as a
%               vector.
% m         =   the number of nodes as considered in the x-direction.
% n         =   the number of nodes as considered in the z-direction.
% Output:
% u_new     =   the solution of the homogenisation problem stored as a
%               matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise matrix form
u_new = zeros(m,n);

% Convert the storage of the solution u
for i = 1:n-1
    u_new(1:m-1,i) = u_old(1+(i-1)*(m-1):i*(m-1))';
end

% Apply periodic boundary condition
u_new(m,:) = u_new(1,:);
u_new(:,n) = u_new(:,1);

end