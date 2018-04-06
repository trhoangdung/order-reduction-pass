function I=inhull(lb,ub);
% INHULL - Creates a interval hull.
% 
% I=inhull(lb,ub)
%
% Creates the interval hull I such that x in I iff x(i) in [lb(i) ub(i)] 
%
% Inputs: 
%   lb - n dimensional column vector
%   ub - n dimensional column vector
%
% Outputs:
%   I  - inhull        
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also INHULL/DOUBLE 

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

I.lb=lb;
I.ub=ub;
I = class(I,'inhull');