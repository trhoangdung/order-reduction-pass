function [M,H] = optimize_trace(P,k)
% This function solve the following problem: 
% Given P = P^T > 0  , find a full row rank matrix H to minize the trace of 

% [I H^T]* [P1 P2; P3 P4]* [I; H]

% input: 1) symetric positive matrix P 
%        2) rank of the matrix H (number of row of H) (k x k is also the size of the abstraction)

% output: 1) the matrix H 

%         2) the corresponding minimized matrix M = P1 + P2*H + P3*H^T +
%         P4*H*H^T

% This optimization problem is the core problem that is used to find an optimal set of initial
% states of the abstraction 
 
% Writer: Dung Tran 9/9/2016

[nP,mP] = size(P); 
n = nP - k; % (n x n) this is the size of original system 
P1 = P(1:n, 1:n);
P2 = P(1:n,n+1:nP);
P3 = P(n+1:nP,1:n); 
P4 = P(n+1:nP,n+1:nP);

I = eye(n); 
h = -inv(eye(n*k)+0.5*kron(I,P4))*vec(P3); 
H = invvec(h,k);

M = P1 + P2*H + P3*H.' +H.'*P4*H; 

end

function X = kron(A,B)
%KRON Kronecker product.
%   kron(A,B) returns the Kronecker product of two matrices A and B, of 
%   dimensions I-by-J and K-by-L respectively. The result is an I*K-by-J*L
%   block matrix in which the (i,j)-th block is defined as A(i,j)*B.

%   Version: 06/02/2011
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)

[I J] = size(A);
[K L] = size(B);

if ~issparse(A) && ~issparse(B)
    
    % Both matrices are dense.
    A = reshape(A,[1 I 1 J]);
    B = reshape(B,[K 1 L 1]);
    X = reshape(bsxfun(@times,A,B),[I*K J*L]);
    
else
    
    % One of the matrices is sparse.
    [ia,ja,sa] = find(A);
    [ib,jb,sb] = find(B);
    ix = bsxfun(@plus,K*(ia(:)-1).',ib(:));
    jx = bsxfun(@plus,L*(ja(:)-1).',jb(:));
    
    % The @and operator is slightly faster for logicals.
    if islogical(sa) && islogical(sb)
        X = sparse(ix,jx,bsxfun(@and,sb(:),sa(:).'),I*K,J*L);
    else
        X = sparse(ix,jx,double(sb(:))*double(sa(:).'),I*K,J*L);
    end
    
end

end

function [x] = vec(X)
% this function implement the vec operator 
x = X(:);
end

function [H] = invvec(h,k)
%This function reconstruct a matrix H from its vec 
% k is the number of row of the matrix
% h is the vector vec 

m = size(h,1)/k; 
H = zeros(k,m);

for i = 1:m
  H(:,i) = h(k*(i-1)+1:k*i);
end

end

