function [X0] = init_set(lb,ub)
% This function follows original simulation approach proposed by ZhiHan


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

end

