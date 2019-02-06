function [A,b1,b2,r] = Construct_A_b(x,z,del,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct_A_b build a linear system Au = b for solving the homogenisation
% problem numerically, where u is the solution to the problem.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% x         =   a vector containing the x coordinates of each nodes.
% z         =   a vector containing the z coordinates of each nodes.
% del       =   a cell containing 2 vectors in the following order:
%               {delx, delz}, where
%               - delx is a vector containing the x-distance between nodes
%               - delz is a vector containing the z-distance between nodes.
% r         =   the type of the material of interest, where
%               - 'a' is the Alluvium Aquifer
%               - 'c' is the Walloon Coal Measure
%               - 't' is the benchmark test case.
% Outputs:
% A         =   A coefficient matrix of the system Ax = b, where x contains
%               the solution to the discretised PDE.
% b1        =   the b vector for j = 1.*
% b2        =   the b vector for j = 2.*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% * 'j' referenced above is not related to j as defined below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialise variables
m = length(x)-1;
n = length(z)-1;
A = zeros(m*n);
b1 = zeros(m*n,1);
b2 = b1;

% Start building process
for i = 1:m
    for j = 1:n
        P = (j-1)*m + i; % Mapping node of reference
        
        % Coordinate mapping for delz
        j_N = j;
        if j == 1
            j_S = n;
        else
            j_S = j-1;
        end
        
        % Coordinate mapping for delx
        i_E = i;
        if i == 1
            i_W = m;
        else
            i_W = i-1;
        end
        
        % Define the weighting terms
        delz_N = del{2}(j_N);
        delx_E = del{1}(i_E);
        delz_S = del{2}(j_S);
        delx_W = del{1}(i_W);
        
        Delz_N = delz_N/2;
        Delx_E = delx_E/2;
        Delz_S = delz_S/2;
        Delx_W = delx_W/2;
        
        % Define variables used for coodinate at the face
        x_P = x(i);
        z_P = z(j);
        x_e = x_P + Delx_E;
        x_w = x_P - Delx_W;
        z_n = z_P + Delz_N;
        z_s = z_P - Delz_S;
        
        switch i
            case 1
                if j == 1 % Lower left corner
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,1+(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K(1+(x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(1+x_w,1+(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = m + 1;
                    E = 2;
                    S = (n-1)*m + 1;
                    W = m;
                elseif j == n % Upper left corner
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K(1+(x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(1+x_w,(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = 1;
                    E = (n-1)*m + 2;
                    S = (n-2)*m + 1;
                    W = n*m;
                else % Left edge
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K(1+(x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(1+x_w,(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = j*m + 1;
                    E = (j-1)*m + 2;
                    S = (j-2)*m + 1;
                    W = j*m;
                end
            case m
                if j == 1 % Lower right corner
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,1+(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K((x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(x_w,1+(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = 2*m;
                    E = 1;
                    S = n*m;
                    W = m - 1;
                elseif j == n % Upper right corner
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K((x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(x_w,(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = m;
                    E = (n-1)*m + 1;
                    S = (n-1)*m;
                    W = n*m - 1;
                else % Right edge
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K((x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(x_w,(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = (j+1)*m;
                    E = (j-1)*m + 1;
                    S = (j-1)*m ;
                    W = j*m - 1;
                end
            otherwise
                if j == 1 % Lower edge
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,1+(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K((x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(x_w,1+(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = m + i;
                    E = i+1;
                    S = (n-1)*m + i;
                    W = i-1;
                elseif j == n % Upper edge
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K((x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(x_w,(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = i;
                    E = (n-1)*m + i+1;
                    S = (n-2)*m + i;
                    W = (n-1)*m + i-1;
                else % Internal node
                    
                    % Evaluate K at each face
                    K_f11 = K(x_e,(z_P + z_s)/2,r);
                    K_f12 = K(x_e,(z_P + z_n)/2,r);
                    K_f21 = K_f12;
                    K_f22 = K((x_P + x_w)/2,z_n,r);
                    K_f31 = K_f22;
                    K_f32 = K(x_w,(z_P + z_s)/2,r);
                    K_f41 = K_f32;
                    K_f42 = K_f11;
                    
                    % Nodes coordinate re-mapping
                    N = j*m + i;
                    E = (j-1)*m + i+1;
                    S = (j-2)*m + i;
                    W = (j-1)*m + i-1;
                end
        end
        
        % Average K at the face
        K_N = Delx_E*K_f22 + Delx_W*K_f21;
        K_E = Delz_N*K_f11 + Delz_S*K_f12;
        K_S = Delx_E*K_f41 + Delx_W*K_f42;
        K_W = Delz_N*K_f32 + Delz_S*K_f31;
        
        % Generate Pth row of A and Pth element of b1 and b2
        A(P,N) = K_N/delz_N;            % North node
        A(P,E) = K_E/delx_E;            % East node
        A(P,S) = K_S/delz_S;            % South node
        A(P,W) = K_W/delx_W;            % West node
        A(P,P) = -sum(A(P,:));          % P node
        
        b1(P) = (K_W - K_E);
        b2(P) = (K_S - K_N);
    end
end

end