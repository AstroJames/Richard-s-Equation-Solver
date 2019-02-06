function [x,residual, k] = precondition_gmres(x0, A, b, tol, precondition,maxiters,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% precondition_gmres performs the generalised miniminal residual method on a linear system.
% The code allows for the input of a preconditioned matrice. If no prconditioning is required the code will
% utilise the identity matrix as the precondition matrice. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% x0           = initial guess vector of x the same size as b
% A            = symmetric nxn matrix for the system Ax=b
% b            = nx1 vector 
% tol          = tolerance level 
% precondition = precondition matrix type. Input either 'no_pre', 'jacobi',
%               'gauss-siedel, 'SOR' or 'incomplete'
% maxiters     = the maximum number of iterations
% r            = material type 'a', 'c' or 't'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
% x            = the nx1 solution vector to Ax = b
% residual     = the residual of the iterative solution
% k            = the number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% determine optimal omega which can be utilised for SOR
if r == 'a' 
    omega = 0.079114; % alluvium omega
    
else 
    omega = 0.4587; % walloon coal omega
end


% Setting up precondition matrix
switch precondition
    case 'no_pre'
        M = eye(length(A(:,1)), length(A(:,1)));
    case 'jacobi'
        M = diag(diag(A));
    case 'gauss-siedel'
        M = tril(A);
    case 'SOR'
        M = diag(diag(A)/omega) + tril(A,-1); 
    case 'incomplete'
        M = ilu(A); % FIX 
    otherwise
        ('You have not selected a correct input. Please choose either no_pre,jacobi, gauss-sidel, SOR or incomplete')
end

%intialise variables 
x = x0;
r0 = b -A*x;
beta = norm(r0);
v1 = r0/beta;
V = zeros(size(A,1),size(A,1));
H = zeros(size(A,1)+1,size(A,1));
V(:,1) = v1;
residual = 20;
tolerance = tol *norm(r0);
k = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while residual > tolerance && k <= maxiters
    
    k = k +1; %iterations counter
    
    % Perform Arnoldi Algorithim to generate AM^(-1)Vm = V(m+1)Hm
    
    u =  M\V(:, k);
    
    V(:, k+1) = A*u;
    
    for j = 1:k
        H(j,k)   = V(:,j)'*V(:,k+1) ;
        V(:,k+1)= V(:,k+1) - H(j,k)*V(:,j);
    end
    
    H(k+1,k) = norm(V(:,k+1));
    
    if H(k+1,k) ~= 0
        V(:,k+1) = V(:,k+1)/ H(k+1, k);
    else
        break
    end
    
    e1  = [1; zeros(k,1)]; % create basis vector

    y   = H(1:k+1, 1:k)\(beta*e1);     % solve linear system for ym  
   
    residual = norm(beta*e1 - H(1:k+1,1:k)*y);  %compute residual 
    
end

u = M\V(:,1:k)*y; %solve linear system for u 


x = x0 + u; % solution of x is linear combination of x0 and u
 
end % end function





