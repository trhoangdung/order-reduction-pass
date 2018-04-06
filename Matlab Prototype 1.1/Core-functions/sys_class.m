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
    
    methods
    end
    
end

