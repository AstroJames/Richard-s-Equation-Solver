function X = TDMA_solve(A,B)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves a TDMA System AX=B using LU factorization
%
% INPUTS:
% A = The matrix of the system Ax = b
% B = The matrix of b vetors for the system Ax = b
% OUTPUTS:
% X = Matrix of Soltions to the system Ax = b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Variables
n = length(A{2});
X=B;
A = TDMA_factor(A);

%U system
for i = 2:n
    X(i,:)=X(i,:)-A{3}(i-1)*X(i-1,:);
end


%L system
X(n,:)=X(n,:)/A{2}(n);
for i = n-1:-1:1
    X(i,:)=(X(i,:)-A{1}(i)*X(i+1,:))/A{2}(i);
end

end



function A = TDMA_factor(A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Factors a tridiagonal matrix A into L and U
%
% INPUTS:
% A = The Matrix that is to be factor
% OUTPUTS:
% A = Matrix of combined L and U, L is in the lower half and U is in the
% Upper half
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=length(A{2});


A{3}(1)=A{3}(1)/A{2}(1); %Starting Lower
for j = 2:n-1
    A{2}(j)=A{2}(j)-A{3}(j-1)*A{1}(j-1); %Diagonal
    A{3}(j)=A{3}(j)/A{2}(j); %Lower
end

A{2}(n)=A{2}(n)-A{3}(n-1)*A{1}(n-1); %Last Diagonal

end
