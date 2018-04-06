function [Rx,Ru]=reach_set(S,dt,N)
% REACH_SET - Computes the reachable set of a CLS
%
% [Rx,Ru]=reach_set(S,dt,N)
%
% Returns two arrays of zonotopes encoding the reachable set 
% of the CLS S computed using a time step dt and N iterations.
% Use PLOT_REACH to display.
% The current version does not handle polytope constraints.
%
% Inputs:
%   S      - CLS
%   dt     - positive number
%   N      - integer 
%
% Outputs:
%   Rx     - array of zonotopes
%   Ru     - array of zonotopes
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also PLOT_REACH 

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




[A,B,C,U,I]=double(S);

if isa(U,'polytope')
  error('REACH_SET: current version does not handle polytope constraints!');
  return;
elseif isa(U,'inhull')
  U=zonotope(U);
end

if isa(I,'polytope')
  error('REACH_SET: current version does not handle polytope constraints!');
  return;
elseif isa(I,'inhull')
  I=zonotope(I);
end

sx=size(A);
su=size(B);
M = [[A , B] ; zeros(su(2),sx(2)+su(2))];
eM = expm(dt*M);
eA= eM(1:sx(1),1:sx(2));
eB= eM(1:sx(1),sx(2)+1:sx(2)+su(2));

Rx=projection(I,C);
Ru=[];
Ut=projection(U,eB);
Xt=projection(I,eA);
for i=1:N 
  Rx=[Rx,projection(Xt,C)];
  Ru=[Ru,projection(Ut,C)];
  Ut=projection(Ut,eA);
  Xt=projection(Xt,eA);
end