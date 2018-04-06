function [d]=eval_prec(S,H,method)
% EVAL_PREC - Evaluates the precision of the approximate bisimulation between a 
%             CLS and its projection
%
% d=eval_prec(S,H,method)
%
% Evaluates the precision of the approximate bisimulation between a CLS S 
% and the projection of S using the surjective map x2 = Hx1.
%
% S must be a CLS under the unstable/stable decomposition form: 
%       | dx/dt = A x + B u, u in U, x(0) in I
%   S:  |
%       | y     = C x
% with A = [ Au 0 ; 0 As ] where Au is an unstable matrix and As is a stable matrix 
% and  C = [ Cu Cs ].
% This decomposition can be obtained using the function DECOMPOSITION.
%
% The matrix H must be of the form H = [ Hu 0 ; 0 Hs ] and the such that the 
% projection of the unstable subsytem of S and the unstable subsystem of S are 
% exactly bisimilar: Cu = Cu*pinv(Hu)*Hu and Hu*Au = Hu*Au*pinv(Hu)*Hu. 
% (Hu = I is obviously an admissible choice). In addition, Hs must be
% chosen such that Hs*As*pinv(Hs) is stable.
%
% method determines the method used to compute the bisimulation function:
% - if method='LMI', the bisimulation function is computed by solving a set of LMIs.
%   This method is accurate but can take some time (recommended for small
%   and medium scale systems)
% - if method='Lyap',  the bisimulation function is computed by solving a Lyapunov equation.
%   This method is not so accurate but is much faster (recommended for large scale systems)
%
% Inputs:
%   S      - CLS of dimension n
%   H      - p x n matrix of rank p
%   method - string  which is either 'LMI' or 'Lyap' 
%
% Outputs:
%   d      - number
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also CLS, CLS/PROJECTION, DECOMPOSITION, BISIMF 

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

sh=size(H);
if sh(1)~= rank(H)
  error('EVAL_PREC: H must be a surjective map!');
  return;
end

[A1,B1,C1,U1,I1]=double(S);

[M,alpha]=bisimf(S,H,method);
  
% computation of beta
G=[eye(size(A1)); H];
Q = G'*M*G;
if isa(I1,'polytope')
  [L,k] = double(I1);
elseif isa(I1,'inhull')
  [lb,ub] = double(I1);
  L=[eye(length(lb));-eye(length(lb))];
  k= [ub;-lb];
elseif isa(I1,'zonotope')
  [c,G] = double(I1);
  si=size(G);
  if si(1)==si(2) & rank(G)==si(2)
    Gp=inv(G);  
    L=[Gp;-Gp];
    k=[ones(si(1),1)+Gp*c; ones(si(1),1)-Gp*c];
  else
    R=enc_orh(I1);
    [c,G] = double(R);
    Gp=inv(G);  
    L=[Gp;-Gp];
    k=[ones(si(1),1)+Gp*c; ones(si(1),1)-Gp*c];
  end
end
su=size(A1);
[x,beta2]=quadprog(-Q-Q',zeros(1,su(2)),L,k);
%[x,mult,how,exitflag,beta2]=mpt_solveQP(-Q-Q',zeros(1,su(2)),L,k);
beta=real(sqrt(-beta2));
  
% precision
d=max(alpha,beta);
