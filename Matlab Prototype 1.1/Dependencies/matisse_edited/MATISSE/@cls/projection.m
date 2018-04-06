function S2=projection(S1,H,type);
% PROJECTION - Projection of the CLS S1
%
% S2=projection(S1,H,type)
%
% Returns the CLS S2 obtained from S1 using the surjective map x2 = H x1.  
% If the set of initial states is an interval hull, argument type
% determines wether its projection is a polytope or a zonotope. 
%
%
% Inputs:
%   S1   - CLS of dimension n
%   H    - m x n matrix of rank m
%   type - string which is either 'polytope' or 'zonotope'
%
% Outputs:
%   S2   - CLS
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also CLS, BISIMF, EVAL_PREC 

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

sh=size(H);
if sh(1)~= rank(H)
  error('PROJECTION: H must be a surjective map!');
  return;
end

[A1,B1,C1,U,I1]=double(S1);
Hp=pinv(H);
A2=H*A1*Hp;
B2=H*B1;
C2=C1*Hp;

if isa(I1,'polytope')
  [G,f]=double(I1);
  I2=polytope(G*[Hp null(H)],f);
  I2=projection(I2,[1:sh(1)]);
elseif isa(I1,'inhull')
  I2=projection(I1,H,type);  
elseif isa(I1,'zonotope')
  I2=projection(I1,H);  
end

S2=cls(A2,B2,C2,U,I2);