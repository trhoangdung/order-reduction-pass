function [sys_bisim,err_bisim] = run_bisimulation(sys,k)
%This function is used to find approximate bisimulation abstraction for
%different order k

[mk,nk] = size(k);
U = inhull(sys.u_lb,sys.u_ub);
I = inhull(sys.x0_lb,sys.x0_ub);
S = cls(sys.sys.a,sys.sys.b,sys.sys.c,U,I);

[mA,nA] = size(sys.sys.a); 

if mA <= 100
    
    for i = 1:mk
      T = tic;
      [sys_bisim(i),err_bisim.err(i)]= reduction(S,k(i),'LMI','appbisim','zonotope');
      %[sys_bisim(i),err(i)]= reduction(S,k(i),'Lyap','appbisim','zonotope');
      %[sys_bisim(i),err(i)]= reduction(S,k(i),'LMI','baltrunc','zonotope');
      %[sys_bisim(i),err(i)]= reduction(S,k(i),'Lyap','baltrunc','zonotope');
      ET = toc(T); % execute time measure
      time(i) = ET;
      err_bisim.time(i) = ET;
    end
    
else
    
    for i = 1:mk
       sys_bisim(i) = S; 
       err_bisim.err(i) = 0;
       err_bisim.time(i) = 0;
    end
    
end

end

