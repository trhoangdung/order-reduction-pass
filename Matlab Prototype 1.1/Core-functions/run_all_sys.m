load('all_sys.mat'); 

%MCS system
k = [5;4]; 
% our approach
[mcsr,mcs_err] = run_reduction(mcs_sys,k,'mcs');
% bisimulation approach
[mcs_bisim,mcs_bisim_err] = run_bisimulation(mcs_sys,k);

% helicopter system 
k = [20;16;10]; % order of reduced system
% our approach
[helir,heli_err] = run_reduction(heli_sys,k,'heli'); 
% bisimulation approach
[heli_bisim,heli_bisim_err] = run_bisimulation(heli_sys,k);


% building model
k = [25;15;6];  
% our approach 
[buildr,build_err] = run_reduction(build_sys,k,'build');
% approximate bisimulation approach 
[build_bisim,build_bisim_err] = run_bisimulation(build_sys,k);


% iss system 
k = [25;10];  
% our approach can be applied for iss system - 270 states, 3 inputs, 3 outputs
% approximate bisimulation method can not be applied for iss system 
% simulation-based approach can not be applied for iss system 

[issr,iss_err] = run_reduction(iss_sys,k,'iss');
[iss_bisim,iss_bisim_err] = run_bisimulation(iss_sys,k);


% pde system 
k = [30;20;10;6];   
% our approach can be applied for pde system - 84 states, 1 input, 1 output
% approximate bisimulation method can not be applied for pde system 
% simulation-based approach can not be applied for pde system 
[pder,pde_err] = run_reduction(pde_sys,k,'pde');
[pde_bisim,pde_bisim_err] = run_bisimulation(build_sys,k);

% FOM system 
k = [20;15;10];   
% our approach can be applied for fom system - 1006 states, 1 input, 1 output
% approximate bisimulation method can not be applied for pde system 
% simulation-based approach can not be applied for pde system 

[fomr,fom_err] = run_reduction(fom_sys,k,'fom');
[fom_bsim,fom_bisim_err] = run_bisimulation(fom_sys,k);

% Synchronous motor control system SMCS

k = [5]; 
[mcsr_sub1,mcs_sub1_err] = run_reduction(mcs_sub1,k,'mcs_sub1');
[mcsr_sub2,mcs_sub2_err] = run_reduction(mcs_sub2,k,'mcs_sub2');

save('all_data.mat');