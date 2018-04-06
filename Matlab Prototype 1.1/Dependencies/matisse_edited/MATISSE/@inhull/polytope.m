function P=polytope(I);
% POLYTOPE - Convert interval hull into a polytope.
% 
% P=polytope(I)
%
% Inputs: 
%   I - inhull
% 
% Ouputs:
%   P - polytope
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also INHULL, POLYTOPE

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

lb=I.lb;
ub=I.ub;
n=length(lb);
P=polytope([eye(n);-eye(n)],[ub;-lb]);