function [S2,p]=reduction2(S1,d,method_sim,method_red,type)
% REDUCTION - Returns an approximation of a CLS and the
%             precision of the associated approximate bisimulation.
%
% [S2,p]=reduction(S1,d,method_sim,method_proj,type)
%
% Returns a d-dimensional approximation S2 of the CLS S1 and the 
% precision p of the approximate bisimulation between S1 and S2. 
% The argument method_sim determines which method must be used to compute the simulation 
% function (see help file of function BISIMF). 
% The argument method_red determines which method is used to reduce the CLS. 
% - If method_red='baltrunc', the method used for reduction is the classical
%   model reduction technique of balanced truncation.
% - If method_red='appbisim', the method used for reduction is based on the
%   computation of a bisimulation function between the CLS and itself. 
% If the set of initial states is an interval hull, argument type
% determines wether its projection is a polytope or a zonotope. 
%
% Inputs:
%   S1         - CLS
%   d          - integer
%   method_sim - string  which is either 'LMI' or 'Lyap' 
%   method_red - string  which is either 'baltrunc' or 'appbisim'
%   type       - string which is either 'polytope' or 'zonotope'
%
% Outputs:
%   S2         - CLS
%   p          - number
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also CLS, BISIMF, EVAL_PREC 

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



[A,B,C,U,I]=double(S1);
n=length(A);
nu=sum(eig(A)>=0);

if nu > d 
  error('REDUCTION: dimension d must be greater than the dimension of the unstable subsystem of the CLS S!');
  return;
elseif nu==d
  % we only keep the unstable subsystem
  % computation of the precision of the approximate simulation
  [Sd,Pd]=decomposition(S1); 
  p=eval_prec(Sd,eye(nu,n),method_sim);
  % computation of the approximation
  S2=projection(S1,Pd,type); 
  S2=projection(S2,eye(nu,n),type);
  return;    
end

% the stable dynamics is approximated by a dynamics of dimension d-nu
% decomposition of the system
[Sd,Pd]=decomposition(S1);

%extraction of the stable subsystem
[Ad,Bd,Cd,Ud,Id]=double(Sd);
As=Ad(nu+1:n,nu+1:n); 
Bs=Bd(nu+1:n,:);
Cs=Cd(:,nu+1:n);

if strcmp(method_red,'baltrunc')
  %reduction using balanced truncation
  su=size(Bs); sc=size(Cs);
  sys=ss(As,Bs,Cs,zeros(sc(1),su(2)));
  [sysb,g,T] = balreal(sys);

  % computation of the approximation
  S2=projection(S1,[eye(nu,n);[zeros(n-nu,nu),T]]*Pd,type);
  S2=projection(S2,[eye(nu,n);[zeros(d-nu,nu),eye(d-nu,n-nu)]],type); 

  % computation of the precision of the approximate simulation
  Sp = projection(Sd,[eye(nu,n);[zeros(n-nu,nu),T]],'zonotope');
  p = eval_prec(Sp,[eye(nu,n);[zeros(d-nu,nu),eye(d-nu,n-nu)]],method_sim);
  
elseif strcmp(method_red,'appbisim')
  if isa(I,'polytope') | isa(U,'polytope')
    error('REDUCTION: method appbisim does not handle polytope constraints!');
    return;
  end 
  %reduction using simulation function
  %lambda=-min(real(eig(As)))/100; 
  lambda=min(-real(eig(As)))/1.01; 
  if strcmp(method_sim,'LMI')
    Ms = sdpvar(length(As));
    F = [(Ms > Cs'*Cs), ((As'+lambda*eye(size(As)))*Ms+Ms*(As+lambda*eye(size(As))) < 0)];
    solvesdp(F,trace(Ms));
    Ms = double(Ms);
  elseif strcmp(method_sim,'Lyap')
    CC= Cs'*Cs;
    Qs = (As'+lambda*eye(size(As)))*CC+CC*(As+lambda*eye(size(As)));
    [V,D] = eig(Qs);
    D = diag((diag(D)>0).*diag(D));
    Qs = V*D*V';
    Ms = CC+lyap(As'+lambda*eye(size(As)),Qs);
  else
    error('REDUCTION: invalid choice of method!');
    return;
  end
  [Vs,D]=eig(As);
  k=1;
  while k<=length(D)
    if ~isreal(D(k,k))
      W(:,1)=0.5*(Vs(:,k)+Vs(:,k+1));
      W(:,2)=0.5*real(-i*Vs(:,k)+i*Vs(:,k+1));
      Vs(:,k:k+1)=W;
      k=k+2;
    else
      k=k+1;
    end
  end
  Vsp=inv(Vs);
  D=diag(Vs'*Ms*Vs);
  if isa(I,'inhull')
    Is=projection(I,[zeros(n-nu,nu),Vsp]*Pd,'zonotope');
  else
    Is=projection(I,[zeros(n-nu,nu),Vsp]*Pd);  
  end 
  if isa(U,'inhull')
    Us=projection(U,[zeros(n-nu,nu),Vsp]*Pd*B,'zonotope');
  else
    Us=projection(U,[zeros(n-nu,nu),Vsp]*Pd*B);
  end 
  [li,ui]=double(enc_inhull(Is));
  [lu,uu]=double(enc_inhull(Us));
  %m=max(abs([li';ui';lu';uu']));
  m=max(abs([li';ui';lu'/lambda;uu'/lambda]));
  D=real(sqrt(diag(D))*diag(m));
  [D,j]=sort(-diag(D));
  %[D,j]=sort(diag(D));
  Vs=Vs(:,j);
  Vsp=inv(Vs);
  % computation of the approximation
  S2=projection(S1,[eye(nu,n);[zeros(n-nu,nu),Vsp]]*Pd,type);
  S2=projection(S2,[eye(nu,n);[zeros(d-nu,nu),eye(d-nu,n-nu)]],type); 

  % computation of the precision of the approximate simulation
  Sp = projection(Sd,[eye(nu,n);[zeros(n-nu,nu),Vsp]],'zonotope');
  p = eval_prec(Sp,[eye(nu,n);[zeros(d-nu,nu),eye(d-nu,n-nu)]],method_sim);
  
else
  error('REDUCTION: invalid choice of method!');
  return;
end