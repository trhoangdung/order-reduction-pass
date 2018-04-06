function [sys_red,sys_bal,sys_err] = model_reduction(system,k)

sys_ori = system.sys;
lb = system.x0_lb;
ub = system.x0_ub; 
u_lb = system.u_lb;
u_ub = system.u_ub;
X0 = system.X0;

sys_red = sys_class;
sys_bal = sys_class; 
sys_err = sys_class;

[sysb, g, T, Ti] = balreal(sys_ori);
sys_bal.sys = sysb;
sys_bal.T = T;

[mT,nT] = size(T);
[lb_b, ub_b] = Init_Con_Trans(lb,ub,T,mT); % get initial state transformed for balanced system
sys_bal.x0_lb = lb_b;
sys_bal.x0_ub = ub_b;
sys_bal.u_lb = u_lb;
sys_bal.u_ub = u_ub;
Xb0 = T*X0;
sys_bal.X0 = Xb0;

% get reduced system
S_A1   = horzcat(eye(k),zeros(k,mT-k));
S_A2   = vertcat(eye(k),zeros(mT-k,k));
S_B  = horzcat(eye(k),zeros(k,mT-k));

A_r = S_A1*sysb.a*S_A2;
B_r = S_B*sysb.b;
C_r = sysb.c*S_B';  

lb_r = S_B*lb_b; % get initial state transformed for reduced system
ub_r = S_B*ub_b;
sys_red.x0_lb = lb_r;
sys_red.x0_ub = ub_r; 
sys_red.u_lb = u_lb;
sys_red.u_ub = u_ub;
Xr0 = Xb0(1:k,:);
sys_red.X0 = Xr0;

sysr = ss(A_r,B_r,C_r,0);
sys_red.sys = sysr; 


% create the new system for computing error
A_e = [sysb.a zeros(mT,k); zeros(k,mT) sysr.a]; 
B_e = [sysb.b;sysr.b];
C_e = [sysb.c -sysr.c];

sys_err.sys = ss(A_e,B_e,C_e,0);
sys_err.u_lb = u_lb; 
sys_err.u_ub = u_ub; 

lb_e = [lb_b;lb_r];
ub_e = [ub_b;ub_r];

sys_err.x0_lb = lb_e;
sys_err.x0_ub = ub_e;
sys_err.X0 = [Xb0;Xr0];

sys_err.k = k;
sys_err.g = g;

end

function [lb_r, ub_r] = Init_Con_Trans(lb,ub,T,k_order)
%Init_Con_Trans:  This function helps to transform the initial condition from
%full-order model to reduced-order model
% T is balanced transformation matrix
% k_order is the order of the reduced-order model

lb_r = zeros(k_order,1); ub_r = zeros(k_order,1); % upper bound and lower bound of initial conditions for reduced system
for i = 1:k_order
    [v_l_op, lb_r(i)] = linprog(T(i,:),[],[],[],[],lb,ub); % linprog() is a optimization function
    [v_u_op, ub_r(i)] = linprog(-T(i,:),[],[],[],[],lb,ub);
    ub_r(i) = -ub_r(i); 
end

end