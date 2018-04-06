classdef E_class
    % This class is used to create an array of reduced system error with
    % different order
    properties
        order; % the corresponding order 
        err; % this is another class contain the error information
        selected_bound; % this is the smallest bound delta derived from combining different method : simulation, optimization, theory
        
        % total time for derive output abstraction and error bounds
        time_opt; % total implementing time for optimization method
        time_theo; % total implementing time for theoretical bound
        time_mix; % total implementing time for mixing of optimization and simulation method
        time_sim; % total implementing time for simulation method
        num_sim; % total number of simulation in simulation method
        
        
    end
    
    
end

