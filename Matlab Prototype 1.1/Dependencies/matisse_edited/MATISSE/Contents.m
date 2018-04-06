% MATISSE (Metrics for Approximate TransItion Systems Simulation and Equivalence)
% Version 1.0 September-2005
%
% Constrained Linear Systems (CLS):
% - Constructor and property access
%   * cls             - Default constructor for CLS objects
%   * double          - Function used to access internal properties of a CLS
% - Methods
%   * decomposition   - Returns an unstable/stable decomposition of a CLS
%   * eval_prec       - Evaluates the precision of the approximate bisimulation between a 
%                       CLS and its projection
%   * projection      - Projection of a CLS
%   * reach_set       - Computes the reachable set of a CLS 
%   * reduction       - Returns an approximation of a CLS and the precision of the 
%                       associated approximate bisimulation
%   * bisimf          - Returns a bisimulation function between a CLS and its projection
%
% Interval hull (INHULL):
% - Constructor and property access
%   * inhull          - Default constructor for interval hull objects
%   * double          - Function used to access internal properties of an interval hull
% - Methods
%   * polytope        - Converts an interval hull object into a polytope object
%   * projection      - Projection of an interval hull
%   * zonotope        - Converts an interval hull object into a zonotope object
%
% Zonotope (ZONOTOPE):
% - Constructor and property access
%   * zonotope        - Default constructor for zonotope objects
%   * double          - Function used to access internal properties of a zonotope
% - Methods
%   * enc_inhull      - Returns the smallest interval hull enclosing a zonotope
%   * enc_orh         - Returns an oriented rectangular hull enclosing a zonotope 
%   * plot_reach      - Plots the reachable set of a CLS
%   * projection      - Projection of a zonotope 
%
% Demos:
%   * matisse_demo1    - Example of use of MATISSE toolbox
%   * matisse_demo2    - Example of use of MATISSE toolbox
%   * matisse_demo3    - Example of use of MATISSE toolbox
%
% Copyright (C) 2005 Antoine Girard
% Department of Electrical and Systems Engineering 
% University of Pennsylvania
%

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

