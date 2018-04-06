function [sysr,err] = run_reduction(sys,k,filename)
% This function is used to run order reduction for different selections of
% k - order

[mk,nk] = size(k);
T = zeros(mk,1);

for i = 1:mk
err(i) = E_class; 
err(i).order = k(i);
str = sprintf('_%d',k(i));
Fn  = [filename str];
T1 = tic; 
[sysr(i),sysb,syse] = model_reduction(sys,k(i));
ET= toc(T1);
Error = error_eval(syse); % compute the error bound e1 e2 with different methods
err_bound = select_bound(Error); % select the minimum bound delta from combining different methods (simulation, optimization and theory)
[SP,UP] = safe_transform(sys,err_bound);
[mS,nS] = size(SP);
for j=1:nS
    sysr(i).S(j) = SP(j);% transformed safety specification of the reduced system
    sysr(i).U(j) = UP(j);% transformed unsafe specification of the reduced system
end

[mE,nE]=size(Error);
time_opt = 0;
time_theo = 0;
time_sim = 0;
time_mix = 0;
num_sim = 0;
for j=1:nE
  time_opt = time_opt + Error(j).time_opt;
  time_theo = time_theo + Error(j).time_theo;
  time_mix = time_mix + Error(j).time_mix;
  time_sim = time_sim + Error(j).time_sim;
  num_sim = num_sim + Error(j).num_sim;
end

err(i).err = Error;
err(i).selected_bound = err_bound; 
err(i).time_opt = time_opt + ET; 
err(i).time_theo = time_theo + ET; 
err(i).time_mix = time_mix + ET; 
err(i).time_sim = time_sim + ET;
err(i).num_sim = num_sim;

print_error(Error,Fn);

end
end

