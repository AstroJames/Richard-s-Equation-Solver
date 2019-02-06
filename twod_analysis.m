function [u1,u2] = twod_analysis(A,b1,b2,m,n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% twod_analysis solve the linear system Au = b twice using pseudo-inverse
% function, then convert the solutions to a matrix storage. This function
% utilised pinv function for the purpose of analysis, and is not to be used
% in the full model implementation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% A     =   A coefficient matrix of the system Ax = b, where x contains the
%           solution to the discretised PDE.
% b1    =   the b vector for j = 1.
% b2    =   the b vector for j = 2.
% m     =   the number of nodes in the x-direction.
% n     =   the number of nodes in the z-direction.
% Output:
% u1    =   the numerical solution to the 2D PDE problem in x-direction,
%           stored in an m x n matrix.
% u2    =   the numerical solution to the 2D PDE problem in z-direction,
%           stored in an m x n matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
u1 = vect_to_mat(pinv(A)*b1,m,n);
u2 = vect_to_mat(pinv(A)*b2,m,n);
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