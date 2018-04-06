function I2=projection(I1,H,type);
% PROJECTION - Projection of an interval hull
% 
% I2=projection(I1,H,type)
%
% Returns I2 = H I1 as a zonotope or a polytope
%
% Inputs: 
%   I1   - inhull of dimension n
%   H    - m x n matrix
%   type - string which is either 'polytope' or 'zonotope'
% 
% Outputs:
%   I2   - polytope or zonotope
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also INHULL, ZONOTOPE, POLYTOPE

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
if type=='zonotope'
    I2=projection(zonotope(I1),H);
elseif type=='polytope'
    Hp=pinv(H);
    I2=polytope(I1);
    [G,f]=double(I2);
    I2=polytope(G*[Hp null(H)],f);
    I2=projection(I2,[1:sh(1)]);
end
    
    