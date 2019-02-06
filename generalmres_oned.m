function [x,residual, k] = generalmres_oned(x0, A, b, tol, maxiters,prevA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function implements the GMRES solver for a 1D Newton Problem. The
% preconditioner is the jacobian from the previous time iteration. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% x0           = the initial guess
% A            = the new Jacobian
% b            = vector valued function F(x)
% tol          = tolerance level 
% maxiters     = the maximum number of iterations
% prevA        = the old Jacobian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
% x            = the solution vector
% residual     = the residual of the iterative solution
% k            = the number of iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%intialise variables 
x           = x0;
r0          = TDMAMultiply(A,x);
r0          = b -r0;
beta        = norm(r0);
v1          = r0/beta;
V           = zeros(size(A{2},1),size(A{2},1));
H           = zeros(size(A{2},1)+1,size(A{2},1));
V(:,1)      = v1;
residual    = 20;
tolerance   = tol *norm(r0);
k           = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


while residual > tolerance && k <= maxiters
    
    k = k +1; %iteration counter
    
    % Perform Arnoldi Algorithim to generate AM^(-1)Vm = V(m+1)Hm
    u           = TDMA_solve(prevA,V(:,k));
    V(:,k+1)    = TDMAMultiply(A,u);
  
    for j = 1:k
        H(j,k)      = V(:,j)'*V(:,k+1) ;
        V(:,k+1)    = V(:,k+1) - H(j,k)*V(:,j);
    end
    
    H(k+1,k) = norm(V(:,k+1));
    
    if H(k+1,k) ~= 0
        V(:,k+1) = V(:,k+1)/ H(k+1, k);
    else
        break
    end
    
    e1  = [1; zeros(k,1)]; % create basis vector

    y   = H(1:k+1, 1:k)\(beta*e1);     % solve linear system for ym  
   
    residual = norm(beta*e1 - H(1:k+1,1:k)*y);%compute residual   
end

u = TDMA_solve(prevA, V(:,1:k)*y);
x = x0 + u; % solution of x is linear combination of x0 and u
end


function AX = TDMAMultiply(A,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TDMAMultiply performs matrix vector multiplication for a tridiagonal
% matrix stored as 3 vectors in a cell array.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% A            = a tridiagonal matrix A stored as 3 vectors in a cell array
% X            = a Vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Outputs:
% AX           = the matrix vector product
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n       = length(A{2});
AX      = zeros(n,1);
AX(1) = A{2}(1)*X(1) + A{1}(1)*X(2);

for i = 2:n-1
    AX(i) = A{3}(i-1)*X(i-1) +A{2}(i)*X(i) + A{1}(i)*X(i+1);
end

AX(n) = A{3}(n-1)*X(n-1) + A{2}(n)*X(n);
end


