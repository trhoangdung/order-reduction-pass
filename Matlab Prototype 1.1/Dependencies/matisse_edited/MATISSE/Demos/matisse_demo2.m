% MATISSE_DEMO2 - Example of use of MATISSE toolbox
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%

%   This file is part of MATISSE 1.0 (SEP2005)
%   Last update : September 15, 2005
%   Copyright (C) 2005 Antoine Girard
%   Department of Electrical and Systems Engineering 
%   University of Pennsylvania
%   agirard@seas.upenn.edu
% 
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program; if not, write to the Free Software
%   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA. 

clear; 

disp(['                                                                                         ';...
      ' Thank you for using MATISSE 1.0                                                         ';...
      ' Metrics for Approximate TransItion Systems Simulation and Equivalence.                  ';...
      '                                                                                         ';...
      ' In this example, we consider a ten dimensional constrained linear system                ';...
      ' This system has a four dimensional unstable subsystem.                                  ';...
      '                                                                                         ';...
      ' First, we compute the precision of the approximate bisimulation between the             ';...
      ' original system and its unstable subsystem.                                             ';...
      '                                                                                         ';...
      ' Then, we compute a six dimensional approximation and evaluate the precision             ';...
      ' of the approximate bisimulation between the original system and this approximation.     ';...
      '                                                                                         ';...
      ' Finally, we compute the reachable set of each system on the time interval [0,2].        ';...
      ' We display the reachable sets, and a set modeling the set of unsafe states.             ';...
      ' For the approximations, this set is bloated by the precision of the approximate         ';...
      ' bisimulation: if an approximation is safe then the original system is safe.             ';...
      '                                                                                         ';...
      ' The approximation of dimension 6 allows to conclude that the                            ';...
      ' original system is safe, whereas the unstable subsystem does not.                       ';...
      '                                                                                         ']);
 
input(' Press the Return key to begin ','s');
echo on;
 
 % Definition of the Constrained Linear System S1
 
 % Matrices of the system
 A1 = [ 1     0    -0.4   2     0.24  1.6  -0.6   0     0.54  0    ;...
        0     0.8  -2    -0.3   4    -0.5   0     0.3   0    -0.18 ;...
        0     0     0     4     0     0     0     0     0     0    ;...
        0     0    -4     0.2   0     0     0     0     0     0    ;...
        0     0     0     0    -0.2  -8     0     0     0     0    ;...
        0     0     0     0     8    -0.2   0     0     0     0    ;...
        0     0     0     0     0     0    -0.3   0     0     0    ;... 
        0     0     0     0     0     0     0    -0.7   0     0    ;...
        0     0     0     0     0     0     0     0    -0.8   0    ;...
        0     0     0     0     0     0     0     0     0    -1   ];
 
 B1 = [ 1     0     0     0     0     0     0     0     0     0    ;...
        0     1     0     0     0     0     0     0     0     0    ;...
        0     0     1     0     0     0     0     0     0     0    ;...
        0     0     0     1     0     0     0     0     0     0    ;...
        0     0     0     0     1     0     0     0     0     0    ;...
        0     0     0     0     0     1     0     0     0     0    ;...
        0     0     0     0     0     0     1     0     0     0    ;...
        0     0     0     0     0     0     0     1     0     0    ;...
        0     0     0     0     0     0     0     0     1     0    ;...
        0     0     0     0     0     0     0     0     0     1   ];
  
 C1 = [  1    0    0    0    0    0    0    0    0    0   ;...
         0    1    0    0    0    0    0    0    0    0   ];
 
 % Definition of the constraints on the input
 U1 = inhull([-0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01 -0.01]',...
             [ 0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01]');
 
 % Definition of the constraints on the initial states of the system         
 I1 = inhull([ 0.29 -0.01  0.19  0.19  0.19  0.19  0.19  0.24  0.19  0.19]',...
             [ 0.31  0     0.2   0.2   0.2   0.2   0.2   0.26  0.21  0.21 ]');

 % Definition of the system        
 S1 = cls(A1,B1,C1,U1,I1);     
  
 echo off;
 input(' Press the Return key to continue ','s');
 echo on;
  
 % Computation of a four dimensional approximation S2
 % [S2,p2]=reduction(S1,4,'LMI','appbisim','zonotope');
 
 echo off;
 input(' Press the Return key to start reduction ','s');
 
 disp('           ')
 [S2,p2]=reduction(S1,4,'LMI','appbisim','zonotope');
 
 text=sprintf(' The precision of the approximate bisimulation between S1 and S2 is %0.5g.',p2);
 disp('           ')
 disp(text)

 echo on;
 
 % Computation of a six dimensional approximation S3
 % [S3,p3]=reduction(S1,6,'LMI','appbisim','zonotope');
 
 echo off;
input(' Press the Return key to start reduction ','s');

disp('           ')
[S3,p3]=reduction(S1,6,'LMI','appbisim','zonotope');
 
text=sprintf(' The precision of the approximate bisimulation between S1 and S3 is %0.5g.',p3);
disp('  ')
disp(text)
disp('  ')

input(' Press the Return key to continue ','s');
echo on;
 
 % Computation of the reachable set of S1 on the time interval [0,4] 
 [Rx1,Ru1]=reach_set(S1,0.01,200);
 
 % Computation of the reachable set of S2 on the time interval [0,4] 
 [Rx2,Ru2]=reach_set(S2,0.01,200);
 
 % Computation of the reachable set of S3 on the time interval [0,4] 
 [Rx3,Ru3]=reach_set(S3,0.01,200);
 
 echo off;
 input(' Press the Return key to display the results ','s');


% Display the unsafe sets
Ux1=0.8*ones(1,101)+0.3*cos([0:2*3.147/100:2*3.147]); Uy1=-1*ones(1,101)+0.3*sin([0:2*3.147/100:2*3.147]);
subplot(2,2,1), fill(Ux1,Uy1,'r');

Ux2=0.8*ones(1,101)+(0.3+p2)*cos([0:2*3.147/100:2*3.147]); Uy2=-1*ones(1,101)+(0.3+p2)*sin([0:2*3.147/100:2*3.147]);
subplot(2,2,2), fill(Ux2,Uy2,'r');

Ux3=0.8*ones(1,101)+(0.3+p3)*cos([0:2*3.147/100:2*3.147]); Uy3=-1*ones(1,101)+(0.3+p3)*sin([0:2*3.147/100:2*3.147]);
subplot(2,2,3), fill(Ux3,Uy3,'r');

% Display the results
subplot(2,2,1), plot_reach(Rx1,Ru1,eye(2),'b'), title('Reachable set of the original 10-d system S1'), axis([0.2 2.2 -1.5 0.1]);
subplot(2,2,2), plot_reach(Rx2,Ru2,eye(2),'b'), title('Reachable set of the 4-d approximation S2'), axis([0.2 2.2 -1.5 0.1]);
subplot(2,2,3), plot_reach(Rx3,Ru3,eye(2),'b'), title('Reachable set of the 6-d approximation S3'), axis([0.2 2.2 -1.5 0.1]);

