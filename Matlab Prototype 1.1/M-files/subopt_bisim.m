function [P,abs] = subopt_bisim( sys, k )
% This function produce the sub-optimal abstraction matrices from 
% a stable autonomous high-dimensional system
% Input: 1) sys: the original system 
%        2) k : the order of the abstraction that we want to obtain

% Output: 1) the minimized trace symetric P of the Lyapunov function 
% V = x^T P x
%         2) The abstraction 

% Dung Tran - 9/15/2016

A = sys.A; 
B = sys.B; 
C = sys.C; 



end

 