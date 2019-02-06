function [Kxx,Kzz] = Tensor_solver(u1,u2,x,z,del,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tensor_solver solve the effective conductivities Kxx and Kzz
% simultaneously using bilinear interpolation and numerical integration.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% u1        =   m x n matrix containing numerical solution of u1.
% u2        =   m x n matrix containing numerical solution of u2.
% x         =   m x 1 vector of all the x-coordinates of the grid.
% z         =   n x 1 vector of all the z-coordinates of the grid.
% del       =   a cell containing 2 vectors in the following order:
%               {delx, delz}, where
%               - delx is a vector containing the x-distance between nodes
%               - delz is a vector containing the z-distance between nodes.
% r         =   the type of the material of interest, where
%               - 'a' is the Alluvium Aquifer,
%               - 'c' is the Walloon Coal Measure,
%               - 't' is the benchmark test case.
% Outputs:
% Kxx       =   effective conductivity in x-direction.
% Kzz       =   effective conductivity in z-direction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise variables
Kxx = 0;
Kzz = 0;
m = length(x);
n = length(z);

% Effective Conductivities Solver
for i = 1:m-1
    for j = 1:n-1
        
        % Set up Bilinear Interpolation
        A = zeros(4,4);
        A(:,1) = ones(4,1);
        A(:,2:3) = [x(i),z(j);x(i+1),z(j);x(i),z(j+1);x(i+1),z(j+1)];
        A(:,4) = A(:,2).*A(:,3);
        b1 = [u1(i,j);u1(i+1,j);u1(i,j+1);u1(i+1,j+1)];
        b2 = [u2(i,j);u2(i+1,j);u2(i,j+1);u2(i+1,j+1)];
        
        % Solve Bilinear Interpolant constant
        C1 = A\b1;
        C2 = A\b2;
        
        % Centre node
        x_bar = (x(i+1)+x(i))/2;
        z_bar = (z(j+1)+z(j))/2;
        
        % Calculate/Estimate Parameters
        u1x = C1(2) + C1(4)*z_bar;
        u2z = C2(3) + C2(4)*x_bar;
        Area = del{1}(i)*del{2}(j);
        K_local = K(x_bar,z_bar,r);
        
        % Accumulating Efective Conductivity terms
        Kxx = Kxx + K_local*(u1x+1)*Area;
        Kzz = Kzz + K_local*(u2z+1)*Area;
    end
end

end