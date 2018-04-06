function [E] = error_eval(err_sys)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%[sysr,sysb,syse] = model_reduction(pde_sys,30);
%err_sys = syse;


A_e = err_sys.sys.a;
B_e = err_sys.sys.b;
C_e = err_sys.sys.c;
lb_e = err_sys.x0_lb;
ub_e = err_sys.x0_ub;
u_lb = err_sys.u_lb; 
u_ub = err_sys.u_ub;
gamma = err_sys.gamma;
k = err_sys.k; 
g = err_sys.g;

[mAe,nAe] = size(A_e);
[mBe,nBe] = size(B_e);
[mCe,nCe] = size(C_e); 

for i = 1:mCe
    Err(i) = error_class;
end

% computing theoretical bound 

e1_time_theo1 = tic; % calculate computing time for theoretical bound of e1
% calculate sup |x0|_2 
p1 = zeros(mAe,1);
for i=1:mAe
p1(i) = max([abs(lb_e(i)) abs(ub_e(i))]);    
end
normp1 = norm(p1,2);
E_e1_time_theo1 = toc(e1_time_theo1); 

for i = 1:mCe
    e1_time_theo2 = tic; 
    C1 = C_e(i,:); 
    % theory bound of e1 
    d = eigs(C1'*C1);
    Err(i).e1_theo = sqrt(d(1))*normp1;
    E_e1_time_theo2 = toc(e1_time_theo2); 
    Err(i).e1_time_theo = E_e1_time_theo1/mCe + E_e1_time_theo2; % average computing time of theoretical bound of e1 of the corresponding output i 
end

e2_time_theo = tic; % calculate computing time for theoretical bound of e2
q = 0; 
for i = k+1:mAe-k
    q = q+(2*i-1)*g(i);
end
% calculate norm u 
umax = zeros(nBe,1);
for i = 1:nBe
   umax(i) = max([abs(u_lb(i)) abs(u_ub(i))]);
end
e2_theo = 4*norm(umax,inf)*q; 
E_e2_time_theo = toc(e2_time_theo); % computing time of theoretical bound of e2
for i = 1:mCe
    % theory bound of e2
    Err(i).e2_theo = e2_theo;
    Err(i).e2_time_theo = E_e2_time_theo/mCe; % average computing time for the theoretical bound of e2 of the output i  
end

% calculate the otimization-based bound

if mAe <= 150 % only solve optimization problem for a error system under 100-dimensions
    for i = 1:mCe
    e1_time_opt = tic; % compute the computing time of optimization-based method 
    C1 = C_e(i,:);   
    % bound of e1 using optimization 
    P = sdpvar(mAe,mAe);
    F = [P>0,A_e'*P+P*A_e < 0,C1'*C1 <=P];
    %optimize(F,trace(P));
    solvesdp(F,trace(P)); % old version of yalmip need to use solvesdp
    P1 = (1/2)*double(P); % old version of yalmip, need to use function double
    %P1 = (1/2)*value(P);
    opts = optimoptions('quadprog','Algorithm','active-set','Display','off'); % this for non-convex casse
    [x,fval,eflag,output,lambda] = quadprog(-P1,zeros(mAe,1),zeros(mAe,mAe),zeros(mAe,1),zeros(mAe,mAe),zeros(mAe,1),lb_e,ub_e,lb_e,opts);% X0 = lb_e
    %[x,fval] = quadprog(-P1,zeros(mAe,1),zeros(mAe,mAe),zeros(mAe,1),zeros(mAe,mAe),zeros(mAe,1),lb_e,ub_e); % algorithm 'interior-point-convex' for convex case
    
    e1_opt = sqrt(-fval);
    Err(i).e1_opt = e1_opt;
    E_e1_time_opt = toc(e1_time_opt); 
    Err(i).e1_time_opt = E_e1_time_opt;
    end

else 
    
    for i=1:mCe
        Err(i).e1_opt = 0; % this means the optimization-based method does not work for > 150-dimensons system 
        Err(i).e1_time_opt = 0;
    end
   
end

% compute error bound e1 and e2 using simulation 

X0 = err_sys.X0; 
[mX,nX] = size(X0);

Nt = 1000;   % 1000 samples
t_end = 10;  % simulate for 10 seconds
U = zeros(nBe, Nt);      % a zero input on all inputs 
T = linspace(0, t_end, Nt);
ymax = zeros(nBe,nX);


for i = 1:mCe % number of output
    Err(i).num_sim_e1 = log2(nX)/mCe; % average number of simulation for e1 bound for the correponding output i
    Err(i).num_sim_e2 = 1;          % number of simulation for e2 bound for the corresponding output i
end


% compute bound of e1 using simulation 
if (sum(X0 ~=0) >= 1)    
    e1_time_sim = tic;
    for i=1:nX % nX is number of vertices in the inital set of states
        [Ysim,tsim,Xsim] = lsim(err_sys.sys, U, T, X0(:,i)');
        S = lsiminfo(Ysim,tsim);       
        for j=1:mCe
            ymax(j,i) = max([abs(S(j).Min) abs(S(j).Max)]);
            
        end       
    end
    E_e1_time_sim = toc(e1_time_sim);
    
    for j=1:mCe
        Err(j).e1_sim = (1+gamma)*max(ymax(j,:));
       % Err(j).e1_sim = max(ymax(j,:));
        Err(j).e1_time_sim = E_e1_time_sim/mCe; % average simulation time for e1 bound for output j
    end 
    
else    
    for j=1:mCe
        Err(j).e1_sim = 0;
        Err(j).e1_time_sim = 0; 
    end   
end

% compute bound of e2 using simulation 

for i=1:mCe    
  e2_time_sim = tic;   
  C1 = C_e(i,:); 
  sys1 = ss(A_e,B_e,C1,0); 
  % simulation bound of e2    
  s2 = stepinfo(sys1); 
        
  e2_sim = 0;
  for j = 1:nBe
  e2_sim = e2_sim + s2(j).Peak*umax(j);
  end
  E_e2_time_sim = toc(e2_time_sim);
  
  Err(i).e2_sim = (1+gamma)*e2_sim; 
 % Err(i).e2_sim = e2_sim; 
  Err(i).e2_time_sim = E_e2_time_sim;  
end


 E = Err;
