classdef error_class
    properties
        e1_theo; % theory bound of e1
        e1_sim;  % simulation bound of e1 
        e1_opt; % bound of e1 using optimization method
        e2_theo; % theory bound of e2 
        e2_sim;  % simulation bound of e2
        e2_opt = 0; % bound of e2 using optimization method, not used now          
        
        e1_time_theo; % time for computing theoretical bound of e1
        e1_time_opt; % time for computing e1_opt 
        e1_time_sim; % time for determine bound of e1 using simulation
        
        e2_time_theo; % time for computing theoretical bound of e2
        e2_time_opt = 0; % time for computing optimized bound of e2: not used now
        e2_time_sim; % time for determine bound of e2 using simulation
        
        num_sim_e1; % number of vertices of initial set (= number of simulation need to be done for computing e1's bound)
        num_sim_e2; % number of simulation need to be done for determine e2's bound
           
    end
     
    methods
        
      function r = e_total_theo(obj) % total theory error bound
         r = obj.e1_theo + obj.e2_theo;  
      end
      
      function r = e_total_opt(obj)  % total optimized error bound 
         r = obj.e1_opt + obj.e2_opt;  
      end
      
      function r = e_total_sim(obj)  % total simulation error bound 
         r = obj.e1_sim + obj.e2_sim;  
      end
      
      function r = e_total_mix(obj)  % total of optimized bound of e1 and simulation bound e2
         r = obj.e1_opt + obj.e2_sim;  
      end
      
        
      function r = time_theo(obj)  % time for computing theoretical bounds of e1 and e2 of "each" pair of output
         r = obj.e1_time_theo + obj.e2_time_theo;
      end
      
      function r = time_opt(obj) % time for computing e1 and e2 bound using optimization method
         r = obj.e1_time_opt + obj.e2_time_opt;  
      end
      
      function r = time_sim(obj) % time for computing e1 and e2 bound using simulation method
         r = obj.e1_time_sim + obj.e2_time_sim;  
      end
      
      function r = time_mix(obj) % time for computing e1 by optimization and e2 by simulation 
         r = obj.e1_time_opt + obj.e2_time_sim;  
      end
      
      function r = num_sim(obj) % number of simulation need to be done to determine bounds of e1 and e2 
         r = obj.num_sim_e1 + obj.num_sim_e2;  
      end
      
   
    end
    
end

