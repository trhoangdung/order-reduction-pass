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

