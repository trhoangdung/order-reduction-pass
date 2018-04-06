% MATISSE_DEMO3 - Example of use of MATISSE toolbox
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

clear; %echo on;
disp([  '                                                                               ';...
        ' Thank you for using MATISSE 1.0                                               ';...
        ' Metrics for Approximate TransItion Systems Simulation and Equivalence.        ';...
        '                                                                               ';...
        ' In this example, we consider the heat diffusion equation of a dimensional rod:';...
        '                                                                               ';...
        '   | d T(x,t)/dt = alpha d^2 T(x,t)/dx^2 + u(x,t)                              ';...
        '   | T(0,t) = T(1,t) = 0                                                       ';...
        '   | T(x,0) = 0                                                                ';...
        '                                                                               ';...
        ' T(x,t) is the temperature on the rod at time t, u(x,t) is the heat source,    ';...
        ' we assume that it is concentrate in x=1/3 and constrained in [1,1.1].         ';... 
        ' The output of the system is the temperature at the point x=2/3.               ';...
        '                                                                               ';...
        ' The partial differential equation is discretized in space (101 nodes). This   ';...
        ' 101 dimensional constrained linear system is our original system.             ';...
        '                                                                               ';...
        ' We compute a ten and twenty dimensional approximations and evaluate the       ';...
        ' the precision of the approximate bisimulation between the original system     ';...
        ' and this approximation.                                                       ';...
        '                                                                               '])
 
input(' Press the Return key to begin ','s');

echo on
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %      Definition of the Constrained Linear System S1                  %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Definition of the conductivity of the rod
 alpha=0.01;
 
 % Parameters for the discretization of the heat equation (101 nodes)
 p=34;
 N=3*p-1;
 
 % The matrices of the system S1 are obtained by spatial discretization of the PDE
 echo off;

A1=diag(-2*ones(N,1));
for k=1:N-1
  A1(k,k+1)=1;
  A1(k+1,k)=1;
end
A1=alpha*(N+1)^2*A1;

B1=zeros(N,1); B1(p)=N;

C1=zeros(1,N); C1(2*p)=1;

 echo on;
  
 % Definition of the constraints on the input
 U1= inhull(1,1.1);
 
 % Definition of the constraints on the initial states of the system
 I1 = inhull(zeros(N,1),0.1*ones(N,1));
 
 % Definition of the 200 dimensional system S1       
 S1 = cls(A1,B1,C1,U1,I1); 
 echo off;
 
 disp('   ')
 input(' Press the Return key to continue ','s');

 echo on;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %      Computation of a 10 dimensional approximation S2                %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
 
 echo off;
input(' Press the Return key to start reduction ','s');
disp('           ')
disp(' [S2,p2]=reduction(S1,10,''Lyap'',''appbisim'',''zonotope'')');
disp('    ')
tic;
[S2,p2]=reduction(S1,10,'Lyap','appbisim','zonotope');
t2=toc;
text1=sprintf(' The time needed to compute S2 and the precision of the approximate');
text2=sprintf(' bisimulation between S1 and S2 is %0.5g.',t2);
disp('           ')
disp(text1)
disp(text2)
text=sprintf(' The precision of the approximate bisimulation between S1 and S2 is %0.5g.',p2);
disp('           ')
disp(text)
disp('           ')

input(' Press the Return key to continue ','s');
echo on;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %       Computation of a 20 dimensional approximation S3               %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 echo off;

 input(' Press the Return key to start reduction ','s');
 disp('           ')
 disp(' [S3,p3]=reduction(S1,20,''Lyap'',''appbisim'',''zonotope'')');
 tic;
 [S3,p3]=reduction(S1,20,'Lyap','appbisim','zonotope');
 t3=toc;
 text1=sprintf(' The time needed to compute S3 and the precision of the approximate');
 text2=sprintf(' bisimulation between S1 and S3 is %0.5g.',t3);
 disp('           ')
 disp(text1)
 disp(text2)
 text=sprintf(' The precision of the approximate bisimulation between S1 and S3 is %0.5g.',p3);
 disp('           ')
 disp(text)
 disp('           ')

 input(' Press the Return key to continue ','s');
 echo on;
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %       Computation of the reachable sets                              %
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 echo off;

 disp('    ')
 disp(' [Rx1,Ru1]=reach_set(S1,0.25,400);');
 tic;
 [Rx1,Ru1]=reach_set(S1,0.25,400);
 tr1=toc;
 disp('   ')
 text=sprintf(' The time needed to compute the reachable set of S1 is %0.5g.',tr1);
 disp(text)
 
 disp('    ')
 disp(' [Rx2,Ru2]=reach_set(S2,0.25,400);');
 tic;
 [Rx2,Ru2]=reach_set(S2,0.25,400);
 tr2=toc;
 disp('   ')
 text=sprintf(' The time needed to compute the reachable set of S2 is %0.5g.',tr2);
 disp(text)
 
 disp('    ')
 disp(' [Rx3,Ru3]=reach_set(S3,0.25,400);');
 tic;
 [Rx3,Ru3]=reach_set(S3,0.25,400);
 tr3=toc;
 disp('   ')
 text=sprintf(' The time needed to compute the reachable set of S3 is %0.5g.',tr3);
 disp(text)
 
 
 
 dt=0.25;
 [c,G]=double(Rx1(1));
 Rxt1=zonotope([0.5*dt;c],[[zeros(1,length(G));G] [0.5*dt;0]]);
 Rut1=[];
 
 for k=1:length(Ru1)
     [c,G]=double(Rx1(k+1));
     Rxt1=[Rxt1 zonotope([(k+0.5)*dt;c],[[zeros(1,length(G));G] [0.5*dt;0]])];
     [c,G]=double(Ru1(k));
     Rut1=[Rut1 zonotope([0;c],[zeros(1,length(G));G])];
 end

 [c,G]=double(Rx2(1));
 Rxt2=zonotope([0.5*dt;c],[[zeros(1,length(G));G] [0.5*dt;0]]);
 Rut2=[];
 
 for k=1:length(Ru2)
     [c,G]=double(Rx2(k+1));
     Rxt2=[Rxt2 zonotope([(k+0.5)*dt;c],[[zeros(1,length(G));G] [0.5*dt;0]])];
     [c,G]=double(Ru2(k));
     Rut2=[Rut2 zonotope([0;c],[zeros(1,length(G));G])];
 end
 
 subplot(2,1,1);
 plot_reach(Rxt1,Rut1,eye(2),'b'), title('Evolution of the reachable set of the original system S1'), axis([0 101 0 14]);
 subplot(2,1,2);
 plot_reach(Rxt2,Rut2,eye(2),'r'), title('Evolution of the reachable set of the 10 dimensional approximation S2'), axis([0 101 0 14]);