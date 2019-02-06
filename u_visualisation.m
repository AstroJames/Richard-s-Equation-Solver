function u_visualisation(x,z,u1,u2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u_visulisation produces visualisation for the u1 and u2 solutions across
% the alluvium and walloon coal unit cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% u1        =   m x n matrix containing numerical solution of u1.
% u2        =   m x n matrix containing numerical solution of u2.
% x         =   m x 1 vector of all the x-coordinates of the grid.
% z         =   n x 1 vector of all the z-coordinates of the grid.
% Outputs:
% Two figures showing an interpolation of the solution of u1 and u2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise variable
[X, Z] = meshgrid(x,z) ; % set up mesh
u1 = u1';
u2 = u2';

% Convert into coordinates for triangulisation
k1 = ~isnan(u1);
x = X(k1);
z = Z(k1);
tri = delaunay(x,z); % set up triangulisation
soln1 = u1(k1);
soln2 = u2(~isnan(u2));

% Plotting the visualisation

figure

subplot(1,2,1) % u1 subplot
trisurf(tri, x, z, soln1) % plot for u1
xlabel('x')
ylabel('z')
view(2)
shading interp
colorbar
axis image
grid off

subplot(1,2,2) % u2 subplot
trisurf(tri, x, z, soln2) % plot for u2
xlabel('x')
ylabel('z')
view(2)
shading interp
colorbar
axis image
grid off

end