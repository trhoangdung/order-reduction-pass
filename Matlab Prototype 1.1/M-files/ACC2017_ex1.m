
% building model 

% create building system

load('build.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 48 states, 1 input, 1 output
R0 = 100*eye(48);
Ru = eye(1);

%A_bs = [-0.0075  5.2756 -0.0010 -0.6634; -5.2756 -0.8574 0.0904 0.9217; 
%    -0.0010  -0.0904  -0.0001  -13.5419; 0.6634 0.9217 13.5419 -1.0040];
%B_bs = [0.0061;0.0645;0.0007;-0.0622];
%C_bs = [-0.0061 0.0645 -0.0007 -0.0622];

A1 = [0 0.7;-0.2 -0.6];
B1 = [0.2;-0.4];
C1 = [1 1];
sys1 = ss(A1,B1,C1,0);



A2 = [-0.6 0.4;-0.7 0.2];
B2 = [-0.6;0.4];
C2 = [1 -1];
sys2 = ss(A2,B2,C2,0);


[mA1,nA1] = size(A1);
[mB1,nB1] = size(B1);
[mC1,nC1] = size(C1);

[mA2,nA2] = size(A2);
[mB2,nB2] = size(B2);
[mC2,nC2] = size(C2);

R0_1 = 100*eye(mA1);
Ru_1 = eye(nB1); 

R0_2 = 100*eye(mA2);
Ru_2 = eye(nB2); 


A_aug = [A1 zeros(mA1,nA2);zeros(mA2,nA1) A2];
B_aug = [B1 zeros(mB1,nB2);zeros(mB2,nB1) B2];
C_aug = [C1 -C2]; 

m = mA1+mA2;
P = sdpvar(m);
M1 = sdpvar(m,m,'full');
M2 = sdpvar(m,m,'full');


a = 1; % alpha
%l1 = 0.5;
%l2 = 0.5;
l1 = sdpvar(1,1); % lamda1
l2 = sdpvar(1,1); % lamda2
e  = sdpvar(1,1); % epsilon

R0_aug = [l1*R0_1 zeros(mA1,nA2);zeros(mA2,nA1) (1-l1)*R0_2]; 
Ru_aug = [l2*Ru_1 zeros(nB1,nB2);zeros(nB2,nB1) (1-l2)*Ru_2]; 

Phi_11 = -transpose(A_aug)*transpose(M1) - M1*A_aug + a*P;
Phi_12 = -M1*B_aug;
Phi_13 = P + M1 - transpose(A_aug)*transpose(M2);
Phi_14 = transpose(C_aug);

Phi_22 = -a*Ru_aug; 
Phi_23 = -transpose(B_aug)*transpose(M2);
Phi_24 = zeros(2*mC1,mC1);

Phi_33 = M2 + transpose(M2);
Phi_34 = zeros(m,mC1);

Phi_44 = -eye(mC1);

%Phi = [Phi_11 Phi_12 Phi_13 Phi_14; transpose(Phi_12) Phi_22 Phi_23 Phi_24; transpose(Phi_13) transpose(Phi_23) Phi_33 Phi_34; transpose(Phi_14) transpose(Phi_24) transpose(Phi_34) Phi_44]; 
 
Phi = [A_aug.'*P+P*A_aug + a*P P*B_aug; B_aug.'*P -a*Ru_aug];


Si = [P C_aug.';C_aug e*eye(mC1)];    

F = [l1 > 0, l1 < 1, l2 > 0, l2 < 1, e > 0, Phi <0, P < R0_aug, Si<0];

solvesdp(F,e);
%solvesdp(F,e,sdpsettings('solver','bmibnb'));
