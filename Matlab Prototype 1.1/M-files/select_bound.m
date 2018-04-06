function [err_bound] = select_bound(Error)
% This function is used to select the best error bound delta from different
% methods (simulation, optimization and theory)
[mE,nE] = size(Error);
e1 = zeros(nE,1);
e2 = zeros(nE,1);
err_bound = zeros(nE,1);
for i=1:nE
    if Error(i).e1_opt ~=0
       e1(i) = Error(i).e1_opt;
    elseif Error(i).e1_sim ~=0
       e1(i) = Error(i).e1_sim;
    else
       e1(i) = Error(i).e1_theo; 
    end 
    
    if Error(i).e2_opt ~=0
       e2(i) = Error(i).e2_opt;
    elseif (Error(i).e2_sim ~=0)&&(Error(i).e2_sim < Error(i).e2_theo)
       e2(i) = Error(i).e2_sim;
    else
       e2(i) = Error(i).e2_theo; 
    end 
    
    err_bound(i) = e1(i) + e2(i);   
end

end

