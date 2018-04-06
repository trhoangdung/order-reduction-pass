load('mcs.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 84 states, 2 inputs, 2 outputs
R0 = 100*eye(6);
Ru = eye(2);
a = 1; 



