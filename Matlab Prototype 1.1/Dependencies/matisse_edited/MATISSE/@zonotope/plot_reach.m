function plot_reach(Rx,Ru,P,S);
% PLOT_REACH - Plots the reachable set of a CLS
%
% plot_reach(Rx,Ru,P,S)
%
% Plots the reachable set of a CLS, Rx and Ru are arrays of zonotopes
% computed by REACH_SET. P is a 2 row matrix determining which projection 
% is used for vizualization of the reachable set. S is a string determining the
% options of the plot (see help file of plot.m for details)
%
% Inputs: 
%   Rx - array of zonotopes
%   Ru - array of zonotopes
%   P  - 2 x p matrix
%   S  - string
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also REACH_SET

%   This file is part of MATISSE 1.0 (SEP2005)
%   Last update : September 14, 2005
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

%Z=Rx(1);
%Z=projection(Z,P);
%R=enc_orh(Z);
%[cr,Gr]=double(R);
%R=[cr+Gr(:,1)+Gr(:,2) cr+Gr(:,1)-Gr(:,2) cr-Gr(:,1)-Gr(:,2) cr-Gr(:,1)+Gr(:,2) cr+Gr(:,1)+Gr(:,2)];
%plot(R(1,:),R(2,:),S);

hold on
c_Ut=zeros(2,1);
G_Ut=[];
M_Ut=zeros(2,2);
for i=1:length(Ru)
  [c_Ut1,G_Ut1]=double(projection(Ru(i),P));
  c_Ut=c_Ut+c_Ut1;
  G_Ut=[G_Ut,G_Ut1];
  [c_Xt,G_Xt]=double(projection(Rx(i+1),P));
 
  c_Rt=c_Xt+c_Ut;
  M_Ut=M_Ut+G_Ut1*G_Ut1';
  M_Rt=M_Ut+G_Xt*G_Xt';
  [Q_Rt,D_Rt]=eig(M_Rt);
  P_Rt=Q_Rt'*[G_Xt,G_Ut];
  D_Rt=diag(sum(abs(P_Rt),2));
  G_Rt=Q_Rt*D_Rt;
  Rt=[c_Rt+G_Rt(:,1)+G_Rt(:,2) c_Rt+G_Rt(:,1)-G_Rt(:,2) c_Rt-G_Rt(:,1)-G_Rt(:,2) c_Rt-G_Rt(:,1)+G_Rt(:,2) c_Rt+G_Rt(:,1)+G_Rt(:,2)];
  plot(Rt(1,:),Rt(2,:),S);  
end
