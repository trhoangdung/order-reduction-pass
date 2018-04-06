function [S2,P]=decomposition(S1)
% DECOMPOSITION - Returns an unstable/stable decomposition of the CLS S1
%
% [S2,P]=decomposition(S1)
%
% Returns a CLS S2 and a bijective map P such that S2 is obtained from S1 
% using the change of variable x2=Px1 and S2 is of the form  
%        | dx/dt = A2 x + B2 u, u in U2, x(0) in I2
%   S2:  |
%        | y     = C2 x
% with A2 = [ Au2 0 ; 0 As2 ] where Au2 is an unstable matrix and As2 is a stable matrix.
%
% Inputs:
%   S1   - CLS of dimension n
%   type - string which is either 'polytope' or 'zonotope'
%
% Outputs:
%   S2   - CLS
%   P    - n x n matrix
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

[A1,B1,C1,U,I1]=double(S1);

if all(eig(A1)<0)
   S2=S1;
   P=eye(size(A1));
   return;
end

% decomposition of A1 
[V1,d1]=eig(A1);
d1=diag(d1);
[l,i]=sort(-real(d1));
d1=diag(d1(i));
V1=V1(:,i); 
V1p=inv(V1);

ds= (d1< 0)*d1;
du= d1-ds;
Au= real(V1*du*V1p);
As= real(V1*ds*V1p);

Es=null(Au);
Ss=size(Es);
Eu=null(As);
Su=size(Eu);

A2=[pinv(Eu)*A1*Eu zeros(Su(2),Ss(2));zeros(Ss(2),Su(2)) pinv(Es)*A1*Es];

% computation of the change of basis matrix
[V2,d2]=eig(A2);
d2=diag(d2);
[l,i]=sort(-real(d2));
d2=diag(d2(i));
V2=V2(:,i);
P=real(V2*V1p);
Pp=inv(P);

% decomposition of S1
B2=P*B1;
C2=C1*Pp;

if isa(I1,'polytope')
  [G,f]=double(I1);
  I2=polytope(G*Pp,f);
elseif isa(I1,'inhull')
  I2=projection(I1,P,'zonotope');  
elseif isa(I1,'zonotope')
  I2=projection(I1,P);  
end

S2=cls(A2,B2,C2,U,I2);