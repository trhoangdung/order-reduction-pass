function Z=zonotope(c,G);
% ZONOTOPE - Creates a zonotope of center c and generators G.
% 
% Z=zonotope(c,G)
%
% Creates a zonotope of center c and generators G: z in Z iff z = c + Gx, x(i) in [-1 1] 
%
% Inputs: 
%   c - n dimensional row vector
%   G - n x q matrix
%
% Outputs:
%   Z - zonotope
%                   
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also ZONOTOPE/DOUBLE

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

Z.c=c;
Z.G=G;
Z=class(Z,'zonotope');