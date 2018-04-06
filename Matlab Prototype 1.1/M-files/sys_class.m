classdef sys_class
    
    properties
        sys; % state space model for the system
        x0_lb; % lower bound of inital state
        x0_ub; % uper bound of intial state
        u_lb; % lower bound of control input
        u_ub; % uper bound of control input
        
        % just used for error system
        k; % order of reduced system just used for error system
        gamma = 0.01; % used to make simulation result become sound 
        g; % gramian of balanced system, used to compute theoretical bound of e2
        
        % just used for balanced system
        T; % balanced transformation matrix
        
        X0;% set of initial states that is used for original simulation method
        
        % note that : one safety specification(S) of the original system we compute one new S
        % and new U(unsafe specification) for the reduced system
        S = safe_spec_class; % safety specification of the system: this is an array of safe_spec_class objects
        U = safe_spec_class; % the unsafe specification of the reduced system (only used for the transformed SP of reduced system)
   
    end
    
    %methods
        %function R0 = get_R0(x0_lb,x0_ub) % this function obtains the smallest ellipse containing all initial states
                
        %end
        %function Ru = get_Ru(u_lb,u_ub)  % this function obtains the smallest ellipse containing the control input set
            
       % end    
    %end
    
  methods
       function [V] = get_vertice(lb,ub)
       % This method obtains vertices of the initial set from the lower bound and
       % uper bound vector 
       [ms,ns] = size(lb);
       % calculate number of vertices in the initial set X0
       pos = zeros(ms,1);
       n = 0;
       for i = 1:ms
           if (lb(i)==ub(i))
           pos(i) = 0;
           else 
           pos(i) = i; 
           n = n+1;        
           end
       end

 if n <= 12 % we will do simulation to compute e1 if the number of vertices is smaller than 2^12
         N = 2^n; % number of vertices in the inital set
         T1 = zeros(n,1); 
         k = 0; 
      for i=1:ms
          if(pos(i)~=0)
          k = k+1; 
          T1(k) = pos(i);
          end
      end
      T2 = zeros(N,n); 
      for i=0:N-1
          T2(i+1,:) = de2bi(i,n); 
      end
      T3 = zeros(N,n);
      for i=1:N
          for j = 1:n
             if T2(i,j)== 0
                T3(i,j) = lb(T1(j));
             else
                T3(i,j) = ub(T1(j));
             end
          end
      end
      X0_set = zeros(ms,N);
      for i = 1:N
           for j = 1:ms
              if pos(j) == 0
                 X0_set(j,i)= lb(j);
              else
                 X0_set(j,i) = 0; 
              end
           end
           for k = 1:n
               X0_set(T1(k),i) = T3(i,k);
           end
    
      end

         X0 = X0_set; 
 else 
         X0 = zeros(ms,1);
 end 
 
 V = X0.'; %  vertices are stored row-wise in MPT toolbox

end

end

    
    
end

