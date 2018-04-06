function [X0] = ori_sim(sys)
% This function follows original simulation approach proposed by ZhiHan


lb = sys.x0_lb; 
ub = sys.x0_ub; 
u_lb = sys.u_lb;
u_ub = sys.u_ub; 

T = sysb.T; 

[mT,nT] = size(T);

% calculate number of vertices in the initial set X0
l = zeros(mT,1);
pos = zeros(mT,1); 
for i = 1:mT
    if (lb(i)==ub(i))
        l(i) = 1;
        
    else 
        l(i) = 0; 
        pos(i) = i;
    end
end
n0 = mT-sum(l); 
N = 2^{n0}; % number of vertices in the inital set




% set x0 for each simulation 
X0_set = zeros(mT,N); % set of inital states including N vertices
e1_max = zeros(N,1); % set of error e1 versus vertices 

   for j = 1:mT
      if l(j) == 1
          for k=1:N
            X0_set(j,k) = lb(j);
          end
      end
   end
   

  
   
 X0_sim = X0_set;
end

