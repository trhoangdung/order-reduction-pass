

% create tranmission line system 

load('tline.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 256 states, 2 inputs, 2 outputs

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states  

for i = 1:100
    lb(i) = -0.01;
    ub(i) = 0.01;
end 

for i = 101:200
    lb(i) = 0.002;
    ub(i) = 0.003;
end 

for i = 201:256
    lb(i) = -0.02;
    ub(i) = 0.01;
end 

u_lb = [0.8;-1;0.5;0.9];
u_ub = [1;1;1;1];

tline_sys = sys_class;
tline_sys.sys = sys; 
tline_sys.x0_lb = lb; 
tline_sys.x0_ub = ub;
tline_sys.u_lb = u_lb; 
tline_sys.u_ub = u_ub; 

% create motor control systems MCS

load('mcs.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 84 states, 2 inputs, 2 outputs

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states 

u_lb = [0.16;0.2];
u_ub = [0.3;0.4]; 

lb = [0.002;0;0;0;0.001;0;0;0]; ub = [0.0025;0;0;0;0.0015;0;0;0]; % initial conditions
X0 = init_set(lb,ub);

mcs_sys = sys_class;
mcs_sys.sys = sys; 
mcs_sys.x0_lb = lb; 
mcs_sys.x0_ub = ub;
mcs_sys.u_lb = u_lb; 
mcs_sys.u_ub = u_ub; 
mcs_sys.X0 = X0;


% create helicopter system 

load('helicopter.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 28 states, 6 inputs, 4 outputs

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states  

lb(1) = 0.1; ub(1) = 0.1;
lb(2) = 0.098; ub(2) = 0.11;
lb(3) = 0.098; ub(3) = 0.102;

for i = 4:8
    lb(i) = 0.1;
    ub(i) = 0.1;
end 

for i = 9:28
    lb(i) = 0;
    ub(i) = 0;
end 
X0 = init_set(lb,ub);

u_lb = [-1;-1;-1;-1;-1;-1];
u_ub = [1;1;1;1;1;1];


heli_sys = sys_class;
heli_sys.sys = sys; 
heli_sys.x0_lb = lb; 
heli_sys.x0_ub = ub;
heli_sys.u_lb = u_lb;
heli_sys.u_ub = u_ub;
heli_sys.X0 = X0;

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

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states 

for i = 1:10
    lb(i) = 0.0002;
    ub(i) = 0.00025;
end

for i = 25:25
    lb(i) = -0.0001;
    ub(i) = 0.0001;
end

X0 = init_set(lb,ub);

u_lb = [0.8];
u_ub = [1];

build_sys = sys_class;
build_sys.sys = sys; 
build_sys.x0_lb = lb; 
build_sys.x0_ub = ub;
build_sys.u_lb = u_lb;
build_sys.u_ub = u_ub;
build_sys.X0 = X0;

% create iss system 
load('iss.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 270 states, 3 inputs, 3 outputs

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states

for i = 1:mA
    lb(i) = -0.0001;
    ub(i) = 0.0001;
end 
X0 = init_set(lb,ub);% can not be used due to exceed the maximum size 
% input constraint of iss 
u_lb = [0;0.8;0.9];
u_ub = [0.1;1;1];
% safety specification of iss 
SP1 = safe_spec_class; 
SP1.name = 'safe';
SP1.type = 'polytope';
V = [-0.0005 0.0005 0.0005;-0.0015 -0.00002 0.0005 ; 0.0005 -0.001 0.0005 ;0.001 0.00002 0.0005;0.0008 0.0006 0.0005;-0.0005 0.0005 -0.0005 ;-0.0015 -0.00002 -0.0005 ; 0.0005 -0.001 -0.0005;0.001 0.00002 -0.0005;0.0008 0.0006 -0.0005];

P = polytope(V);
[H,K] = double(P);

SP1.matrix_A = H;
SP1.matrix_B = -K;

iss_sys = sys_class;
iss_sys.sys = sys; 
iss_sys.x0_lb = lb; 
iss_sys.x0_ub = ub;
iss_sys.u_lb = u_lb; 
iss_sys.u_ub = u_ub; 
iss_sys.X0 = X0;
iss_sys.S=SP1;


% create pde system 

load('pde.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 84 states, 1 inputs, 1 outputs

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states  

for i = 1:64
    lb(i) = 0;
    ub(i) = 0;
end 

for i = 65:80
    lb(i) = 0.001;
    ub(i) = 0.0015;
end 

for i = 81:84
    lb(i) = -0.002;
    ub(i) = -0.0015;
end 
X0 = init_set(lb,ub);% can not be used due to exceed the maximum size 

u_lb = [0.5];
u_ub = [1];

pde_sys = sys_class;
pde_sys.sys = sys; 
pde_sys.x0_lb = lb; 
pde_sys.x0_ub = ub;
pde_sys.u_lb = u_lb; 
pde_sys.u_ub = u_ub; 
pde_sys.X0 = X0;

% create FOM system 

load('fom.mat');
A = full(A); % convert sparse matrix to full matrix
B = full(B);
C = full(C); 
[mA,nA] = size(A);
[mB,nB] = size(B);
[mC,nC] = size(C);
D = zeros(mC,nB);
sys = ss(A,B,C,D); % the original model with 1006 states, 1 input, 1 output

lb = zeros(mA,1); ub = zeros(mA,1); % upper bound and lower bound of initial states  

for i = 1:400
    lb(i) = -0.0001;
    ub(i) = 0.0001;
end 

for i = 401:800
    lb(i) = 0;
    ub(i) = 0;
end 

for i = 801:1006
    lb(i) = 0;
    ub(i) = 0;
end 

X0 = init_set(lb,ub);% can not be used due to exceed the maximum size 

u_lb = [-1];
u_ub = [1];
fom_sys = sys_class;
fom_sys.sys = sys; 
fom_sys.x0_lb = lb; 
fom_sys.x0_ub = ub;
fom_sys.u_lb = u_lb;
fom_sys.u_ub = u_ub;
fom_sys.X0 = X0;

% create two sub model of synchronous switching motor system
% motor parameters 
J = 3.2284E-6;
b = 3.5077E-6;
K = 0.0274;
R = 4;
L = 2.75E-6;

% control with integrator (add disturbance compensator)
Aa = [0 1 0 0
      0 -b/J K/J 0
      0 -K/L -R/L 0
      1 0 0 0];
Ba = [0 ; 0 ; 1/L ; 0 ];
Br = [0 ; 0 ; 0; -1];
Ca = [1 0 0 0;0 1 0 0];
%Ca = [1 0 0 0];
%Ca = [0 1 0 0];
Da = 0;

p1 = -100+100i;
p2 = -100-100i;
p3 = -200;
p4 = -300;
Ka = place(Aa,Ba,[p1,p2,p3,p4]);

A_core = Aa-Ba*Ka;
B_core = Br;
C_core = Ca;
D_core = Da; 

A = blkdiag(A_core,A_core);
B = blkdiag(B_core,B_core);
C = [1 0 0 0 0 0 0 0; 1 0 0 0 -1 0 0 0]; % y1 = x1, y2 = x1 - x5 (position error) 
sub1 = ss(A,B,C,0);
sub2 = ss(A,-B,C,0); % inverse direction of two motors

u_lb = [0.16;0.16]; % control input for motor 1
u_ub = [0.2;0.22]; % control input for motor 2

lb1 = [-0.002;0;0;0;-0.001;0;0;0]; ub1 = [0.0025;0;0;0;0.002;0;0;0]; % initial conditions for mode 1
lb2 = [-0.001;0;0;0;-0.002;0;0;0]; ub2 = [0.001;0;0;0;0.003;0;0;0]; % initial conditions for mode 2

mcs_sub1 = sys_class; 
mcs_sub2 = sys_class; 

mcs_sub1.sys = sub1;
mcs_sub2.sys = sub2; 
mcs_sub1.u_lb = u_lb;
mcs_sub1.u_ub = u_ub;
mcs_sub2.u_lb = u_lb;
mcs_sub2.u_ub = u_ub; 
mcs_sub1.x0_lb = lb1;
mcs_sub1.x0_ub = ub1; 
mcs_sub2.x0_lb = lb2;
mcs_sub2.x0_ub = ub2; 
X01 = init_set(lb1,ub1);%
X02 = init_set(lb2,ub2); 
mcs_sub1.X0 = X01;
mcs_sub2.X0 = X02;

V1 = [0.25 0.12; 0.25 0.2;0.4 0.12; 0.4 0.2];
V2 = -V1; 

P1 = polytope(V1);
[E1,c1] = mpt_getInnerEllipsoid(P1); 
P2 = polytope(V2);
[E2,c2] = mpt_getInnerEllipsoid(P2); 

SP1 = safe_spec_class; 
SP1.name = 'unsafe';
SP1.type = 'ellipsoid';
%SP1.matrix_P = [400 0;0 1600];
SP1.matrix_P = E1;
SP1.radius = 1;
SP1.center = c1;
%SP1.center = [0.35;0.175];

SP2 = safe_spec_class; 
SP2.name = 'unsafe';
SP2.type = 'ellipsoid';
%SP2.matrix_P = [400 0;0 1600];
SP2.matrix_P = E2;
SP2.radius = 1;
SP2.center = c2;
%SP2.center = [-0.35;-0.175];

mcs_sub1.S = SP1;
mcs_sub2.S = SP2; 




save('all_sys.mat');
