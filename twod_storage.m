function [A, b1, b2, p_inverse] = twod_storage(A,b1, b2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% twod_storage takes the matrix A and two vectors b1 and b2 and reorder
% them through RCM ordering. The reordered matrix A then get converted into
% a sparse storage. An inverse permutation vector, p_inverse, is given as
% an additional output to restore the order after the system has been fully
% solved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% A     =   A coefficient matrix of the system Ax = b, where x contains the
%           solution to the discretised PDE.
% b1    =   the b vector for j = 1.
% b2    =   the b vector for j = 2.
% Outputs:
% A             =   A matrix in sparse storage with RCM ordering.
% b1            =   b1 vector with RCM ordering.
% b2            =   b2 vector with RCM ordering.
% p_inverse     =   the inverse permutation vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Obtaining permutation vector
p   = symrcm(A);

% Reorder linear system and convert storage to sparse
A   = A(p,p);
A   = sparse(A);
b1  = b1(p);
b2  = b2(p);

% Calculate inverse permuation
n = length(A(:,1));
I = eye(n);
P = I(p,:);
p_inverse = P\[1:n]';


end