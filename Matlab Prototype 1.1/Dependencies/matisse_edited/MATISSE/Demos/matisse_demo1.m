% MATISSE_DEMO1 - Example of use of MATISSE toolbox
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

disp(['                                                                                     ';...   
      ' Thank you for using MATISSE 1.0                                                     ';...    
      ' Metrics for Approximate TransItion Systems Simulation and Equivalence.              ';...
      '                                                                                     ';...
      ' In this example, we define a ten dimensional constrained linear system              ';...
      ' This system is asymptotically stable.                                               ';...
      '                                                                                     ';...
      ' We compute approximations of dimension 5 and 7 and evaluate the precisions of       ';...
      ' the approximate bisimulations between the original system and its approximations.   ';...
      '                                                                                     ';...
      ' We compute the reachable set of each system on the time interval [0,4].             ';...
      ' We display the reachable sets, and a set modeling the set of unsafe states.         ';...
      ' For the approximations, this set is bloated by the precision of the approximate     ';...
      ' bisimulation: if an approximation is safe then the original system is safe.         ';...
      '                                                                                     ';...
      ' The approximation of dimension 7 allows to conclude that the original               ';...
      ' system is safe, whereas the approximation of dimension 5 does not.                  ';...
      '                                                                                     ']);
 
input(' Press the Return key to begin ','s');
echo on;
  
 % Definition of the Constrained Linear System S1
 
 % Matrices of the system
 A1 = [ -0.1  1    1    0    0    0    0    0    0    0   ;...
        -1   -0.1  0    0    0    2    0    0    0    0   ;...
         0    0   -0.3  1    0    0    0    0    0    0   ;...
         0    0    0   -0.1 -8    0    0    0    0    0   ;...
         0    0    0    8   -0.1  0    0    0    0    0   ;...
         0    0    0    0    0   -0.1  1    1    0    0   ;...
         0    0    0    0    0   -0.1 -0.1  0    0    0   ;...
         0    0    0    0    0    0    0   -0.6 -1    1   ;...
         0    0    0    0    0    0    0    1   -0.6  0   ;...
         0    0    0    0    0    0    0    0    0   -0.1 ];
 
 B1 = [  0    0    0    0    0    0    0    0    0    1   ]';
 
 C1 = [  1    0    0    0    0    0    0    0    0    0   ;...
         0    1    0    0    0    0    0    0    0    0   ];
 
 % Definition of the constraints on the input
 U1 = inhull(-0.05,0.05);

 % Definition of the constraints on the initial states of the system
 I1 = inhull([ 9 0 -0.1 -0.1 -2 -0.1 -0.1 -0.1 -0.1 -0.1]',...
             [10 1  0.1  0.1 -1  0.1  0.1  0.1  0.1  0.1]');
 
 % Definition of the system        
 S1 = cls(A1,B1,C1,U1,I1);     
  
 echo off;
input(' Press the Return key to continue ','s');
echo on;
 
 % Computation of a five dimensional approximation S2
 % [S2,p2]=reduction(S1,5,'LMI','appbisim','zonotope');
 
 echo off;
input(' Press the Return key to start reduction ','s');

disp('           ')
[S2,p2]=reduction(S1,5,'LMI','appbisim','zonotope');

text=sprintf(' The precision of the approximate bisimulation between S1 and S2 is %0.5g.',p2);
disp('           ')
disp(text)

echo on;
 
 % Computation of a seven dimensional approximation S3
 % [S3,p3]=reduction(S1,7,'LMI','appbisim','zonotope');
 
 echo off;
input(' Press the Return key to start reduction ','s');

disp('           ')
[S3,p3]=reduction(S1,7,'LMI','appbisim','zonotope');
 
text=sprintf(' The precision of the approximate bisimulation between S1 and S3 is %0.5g.',p3);
disp('           ')
disp(text)
disp('  ')

input(' Press the Return key to continue ','s');
echo on;
 
 % Computation of the reachable set of S1 on the time interval [0,4] 
 [Rx1,Ru1]=reach_set(S1,0.02,200);
 
 % Computation of the reachable set of S2 on the time interval [0,4] 
 [Rx2,Ru2]=reach_set(S2,0.02,200);
 
 % Computation of the reachable set of S3 on the time interval [0,4] 
 [Rx3,Ru3]=reach_set(S3,0.02,200);
 
 echo off;
input(' Press the Return key to display the results ','s');

% Display the unsafe sets
Ux1=-9*ones(1,101)+2.5*cos([0:2*3.147/100:2*3.147]); Uy1=-9*ones(1,101)+2.5*sin([0:2*3.147/100:2*3.147]);
subplot(2,2,1), fill(Ux1,Uy1,'r');

Ux2=-9*ones(1,101)+(2.5+p2)*cos([0:2*3.147/100:2*3.147]); Uy2=-9*ones(1,101)+(2.5+p2)*sin([0:2*3.147/100:2*3.147]);
subplot(2,2,2), fill(Ux2,Uy2,'r');

Ux3=-9*ones(1,101)+(2.5+p3)*cos([0:2*3.147/100:2*3.147]); Uy3=-9*ones(1,101)+(2.5+p3)*sin([0:2*3.147/100:2*3.147]);
subplot(2,2,3), fill(Ux3,Uy3,'r');

% Display the results
subplot(2,2,1), plot_reach(Rx1,Ru1,eye(2),'b'), title('Reachable set of the original 10-d system S1');
subplot(2,2,2), plot_reach(Rx2,Ru2,eye(2),'b'), title('Reachable set of the 5-d approximation S2');
subplot(2,2,3), plot_reach(Rx3,Ru3,eye(2),'b'), title('Reachable set of the 7-d approximation S3');
