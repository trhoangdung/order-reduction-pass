function S=cls(A,B,C,U,I);
% CLS - Creates a constrained linear system.
% 
% S=cls(A,B,C,U,I)
%
% Creates the constrained linear system S given by:
%        | dx/dt = A x + B u, u in U, x(0) in I
%    S:  |
%        | y     = C x
%
% Inputs: 
%   A - n x n matrix
%   B - n x m matrix
%   C - p x n matrix
%   U - set of inputs is either a polytope, a interval hull or a
%       zonotope of dimension m
%   I - set of initial states is either a polytope, a interval hull or a
%       zonotope of dimension n
%
% Outputs:
%   S - CLS
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also CLS/DOUBLE, INHULL, POLYTOPE, ZONOTOPE 

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

S.A=A;
S.B=B;
S.C=C;
S.U=U;
S.I=I;
S=class(S,'cls');