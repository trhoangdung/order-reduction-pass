function [M,alpha]=bisimf(S,H,method)
% BISIMF - Returns a bisimulation function between a CLS and its projection 
%
% [M,alpha]=bisimf(S,H,m)
%
% Returns a bisimulation function between the CLS S and the projection of S 
% using the surjective map x2 = Hx1.
%
% S must be a CLS under the unstable/stable decomposition form: 
%       | dx/dt = A x + B u, u in U, x(0) in I
%   S:  |
%       | y     = C x
% with A = [ Au 0 ; 0 As ] where Au is an unstable matrix and As is a stable matrix 
% and  C = [ Cu Cs ].
% This decomposition can be obtained using the function DECOMPOSITION.
%
% The matrix H must be of the form H = [ Hu 0 ; 0 Hs ] and such that the 
% projection of the unstable subsytem of S and the unstable subsystem of S 
% are exactly bisimilar: Cu = Cu*pinv(Hu)*Hu and Hu*Au = Hu*Au*pinv(Hu)*Hu 
% (Hu = I is obviously an admissible choice). In addition, Hs must be
% chosen such that Hs*As*pinv(Hs) is stable.
%
% The bisimulation function is then of the form
%        | + infinity                     if [ Hu 0 ] x1 /= [ I 0 ] x2
% V(x) = |
%        | max( sqrt( x'*M*x ) ,  alpha ) if [ Hu 0 ] x1  = [ I 0 ] x2 
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
%   M      - (n+p) x (n+p) matrix
%   alpha  - number
%
% Author:   Antoine Girard,
%           Department of Electrical and Systems Engineering 
%           University of Pennsylvania
%
% see also CLS, CLS/PROJECTION, DECOMPOSITION, EVAL_PREC 

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
  error('BISIMF: H must be a surjective map!');
  return;
end

[A1,B1,C1,U1,I1]=double(S);
nu=sum(eig(A1)>=0);
n=length(A1);

if nu==0 
  % system is stable: no further restriction on H
  Hp=pinv(H);
  G=[eye(size(A1)); H];
  A=[A1 zeros(size(H')); zeros(size(H)) H*A1*Hp];
  if sum(eig(A)>=0)~=0
      error('BISIMF: matrix H does not preserve the stability of the system!');
      return;
  end
  C=[C1 -C1*Hp];
  %lambda=-min(real(eig(A)))/100;
  lambda=min(-real(eig(A)))/1.01;
  %eig(A+eye(size(A))*lambda)
  P=eye(size(A));
  Bs1=B1;
  As1=A1;
elseif nu==n
  % system is unstable  
  error('BISIMF: cls S must have a stable subsystem!');
  return;
else
  % decomposition of the system
  Au1=A1(1:nu,1:nu); 
  Z1=A1(1:nu,nu+1:n);
  Z2=A1(nu+1:n,1:nu);
  As1=A1(nu+1:n,nu+1:n);
  Cu1=C1(:,1:nu);
  Cs1=C1(:,nu+1:n);
  Bu1=B1(1:nu,:);
  Bs1=B1(nu+1:n,:);
  if any(any(Z1)) | any(any(Z2)) | sum(eig(Au1)<0)~=0 | sum(eig(As1)>=0)~=0
      error('BISIMF: cls S is not under unstable/stable decomposition form, see help file for details');
      return;
  end
  Z1=H(:,nu+1:n);
  ku=0;
  while ku < sh(1) & all(Z1(ku+1,:)==0)
    ku=ku+1;  
  end
  if ku==0
    error('BISIMF: matrix H is not under appropriate form, see help file for details');
    return;
  else
    Hu=H(1:ku,1:nu);
    Hup=pinv(Hu);
    if max(max(abs(Cu1-Cu1*Hup*Hu)))> 1e-14 | max(max(abs(Hu*Au1-Hu*Au1*Hup*Hu)))> 1e-14
      error('BISIMF: matrix H is not under appropriate form, see help file for details');
      return;
    end
    if ku==sh(1)
      % the stable dynamics is completely ingored
      % the system is under an appropriate form
      G=eye(size(As1));
      A=As1;
      C=Cs1;
      %lambda=-min(real(eig(A)))/100;
      lambda=min(-real(eig(A)))/1.01;
      P=[zeros(n-nu,nu) eye(n-nu) zeros(n-nu,sh(1))]; 
    else    
      Z2=H(ku+1:sh(1),1:nu);
      Hs=H(ku+1:sh(1),nu+1:n);
      if any(any(Z2))
        error('BISIMF: matrix H is not under appropriate form, see help file for details');
        return;
      end
      % the system is under an appropriate form
      Hsp=pinv(Hs);
      G=[eye(size(As1)); Hs];
      A=[As1 zeros(size(Hs')); zeros(size(Hs)) Hs*As1*Hsp];
      if sum(eig(A)>=0)~=0
        error('BISIMF: matrix H does not preserve the stability of the stable subsystem!');
        return;
      end
      C=[Cs1 -Cs1*Hsp];
      %lambda=-min(real(eig(A)))/100;
      lambda=min(-real(eig(A)))/1.01;
      % lambda = 0.00001;
      P=[[zeros(n-nu,nu) eye(n-nu) zeros(n-nu,sh(1))]; [zeros(sh(1)-ku,n) zeros(sh(1)-ku,ku) eye(sh(1)-ku)]];
    end
  end
end
  
% computation of M
if strcmp(method,'LMI')
  % LMI method
  M = sdpvar(length(A));
  F = [M > C'*C, ((A'+lambda*eye(size(A)))*M+M*(A+lambda*eye(size(A))) < 0)];
  %solvesdp(F,trace((G*G')*M));
  solvesdp(F,trace(G'*M*G)); 
  %solvesdp(F,trace((G*(Bs1*Bs1'/lambda+eye(size(As1)))*G')*M));
  %solvesdp(F,trace((G*(Bs1*Bs1')*G')*M));
  M = double(M);
  %eig(M - C'*C)
  %eig((A'+lambda*eye(size(A)))*M+M*(A+lambda*eye(size(A))))
  Q = G'*M*G;
  M=P'*M*P;
elseif strcmp(method,'Lyap')
  % solving Lyapunov equations
  CC= C'*C;
  Q = (A'+lambda*eye(size(A)))*CC+CC*(A+lambda*eye(size(A)));
  [V,d] = eig(Q);
  V=real(V);d=real(d);
  d = diag((diag(d)>0).*diag(d));
  Qp = V*d*V';
  M = CC+lyap(A'+lambda*eye(size(A)),Qp);
  Q = G'*M*G;
  M=P'*M*P;
else
  error('BISIMF: invalid choice of method m');
  return;
end
  
% computation of alpha
if isa(U1,'polytope')
  [L,k] = double(U1);
elseif isa(U1,'inhull')
  [lb,ub] = double(U1);
  L=[eye(length(lb));-eye(length(lb))];
  k= [ub;-lb];
elseif isa(U1,'zonotope')
  [c,G] = double(U1);
  su=size(G);
  if su(1)==su(2) & rank(G)==su(2)
    Gp=inv(G);  
    L=[Gp;-Gp];
    k=[ones(su(1),1)+Gp*c; ones(su(1),1)-Gp*c];
  else
    R=enc_orh(U1);
    [c,G] = double(R);
    Gp=inv(G);  
    L=[Gp;-Gp];
    k=[ones(su(1),1)+Gp*c; ones(su(1),1)-Gp*c];
  end
end
su=size(B1);
[u,alpha2]=quadprog(-Bs1'*Q*Bs1-(Bs1'*Q*Bs1)',zeros(1,su(2)),L,k)
%[u,multi,how,exitflag,alpha2]=mpt_solveQP(-Bs1'*Q*Bs1-(Bs1'*Q*Bs1)',zeros(1,su(2)),L,k);
alpha=sqrt(-alpha2)/lambda;
  