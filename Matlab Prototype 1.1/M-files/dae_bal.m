function [ dae_sys_red ] = dae_bal(dae_sys,lf)
% This function obtain the reduced-order differential algerbraic system
% from a high-order dae system. 
%   Detailed explanation: 
%   It implements the Algorithm 3.3 in the paper "Balanced Truncation Model Reduction for Large-Scale Systems in Descriptor Form " 
% input: the original dae system, lf is the order of the reduced proper part
% output: the reduced dae system

A = dae_sys.a; 
B = dae_sys.b;
C = dae_sys.c; 
D = dae_sys.d; 
E = dae_sys.e; 

[m,n] = size(A); 
[re] = rank(E); 

% step 1: compute generalized Schur form
[EE,AA,VT,U] = qz(E,A,'real'); 
V = VT.'; 
Ef = EE(1:re,1:re);
Eu = EE(1:re,re+1:m);
Einf = EE(re+1:m,re+1:m);
Af = AA(1:re,1:re);
Au = AA(1:re,re+1:m);
Ainf = AA(re+1:m,re+1:m);

% step 2: compute the matrices V^TB = [Bu;Binf] and CU = [Cf, Cu]
BB = VT*B;
Bu = BB(1:re,:);
Binf = BB(re+1:m,:); 

CC = C*U;
Cf = CC(:,1:re);
Cu = CC(:,re+1:m);

% step 3: Solve the system of generalized Sylvester equations 
%   Ef*Y - Z*Einf = -Eu, 
%   Af*Y - Z*Ainf = -Au. -> Z = (Af*Y+Au)*Ainf^(-1)


% -> Ef*Y - Af*Y*Ainf^(-1)*Einf = -Eu + Au*Ainf^(-1)*Einf 
% -> Af^-1*Ef*Y - Y*Ainf^-1*Einf = Af^-1*(-Eu + Au*Ainf^(-1)*Einf)

Y = lyap((inv(Af))*Ef,-inv(Ainf)*Einf,-inv(Af)*(-Eu + Au*(inv(Ainf))*Einf));
Z = (Af*Y+Au)*(inv(Ainf));

% step 4: Compute the cholesky factor Rf, Lf, Rinf and Linf of the
% solutions Xpc = Rf*Rf^T, Xpo = Lf*Lf^T, Xic = Rinf*Rinf^T and Xio =
% Linf*Linf^T of the generalizized Lyapunov equations:

% (1) Ef*Xpc*Af^T + Af*Xpc*Ef^T = -(Bu-Z*Binf)*(Bu-Z*Binf)^T, 
% (2) Ef^T*Xpo*Af + Af^T*Xpo*Ef = -Cf^T*Cf,
% (3) Ainf*Xic*Ainf^T -Einf*Xic*Einf^T = Binf*Binf^T, 
% (4) Ainf^T*Xio*Ainf - Einf^T*Xio*Einf = (Cf*Y+Cu)^T*(Cf*Y+Cu) 

Xpc = lyap(Ef,(Bu-Z*Binf)*((Bu-Z*Binf).'),[],Af);
Xpo = lyap(Ef.',(Cf.')*Cf,[],Af.');
Xic = dlyap(Ainf,-Binf*(Binf.'),[],Einf);
Xio = dlyap(Ainf.',-((Cf*Y+Cu).')*(Cf*Y+Cu),[],Einf.');

RfT = chol(Xpc);
Rf = RfT.';
LfT = chol(Xpo);
Lf = LfT.';
RinfT = chol(Xic);
Rinf = RinfT.';
LinfT = chol(Xio);
Linf = LinfT.'; 

% step 5: Compute the skinny singular value decompositions 

% Lf^T*Ef*Rf = [U1, U2]*[Z1 0; 0 Z2][V1, V2]^T, 
% Linf^T*Ainf*Rinf = U3*D3*V3^T,

r = rank(Lf.'*Ef*Rf); % r = nf  : the dimension of original proper part
linf = rank(Linf.'*Ainf*Rinf);

[UU,ZZ,VV] = svd(Lf.'*Ef*Rf); 
[U3,D3,V3] = svd(Linf.'*Ainf*Rinf);

% What is lf parameter? lf should be an input of the function 

Z1 = ZZ(1:lf,1:lf);
Z2 = ZZ(lf+1:r,lf+1:r); 

U1 = UU(:,1:lf);
U2 = UU(:,lf+1:r);
V1 = VV(:,1:lf);
V2 = VV(:,lf+1:r);

% Step 6: compute Wf = Lf*U1*Z1^(-1/2), Winf = Linf*U3*D3^(-1/2)
% Tf = Rf*V1*Z1^(-1/2) and Tinf = Rinf*V3*D3^(-1/2)

Wf = Lf*U1*Z1^(-1/2);
Winf = Linf*U3*D3^(-1/2); 
Tf = Rf*V1*Z1^(-1/2);
Tinf = Rinf*V3*D3^(-1/2);

% Step 7: Compute the reduced-order system 

E_r = [eye(lf) zeros(lf,linf); zeros(linf,lf) (Winf.')*Einf*Tinf];
A_r = [(Wf.')*Af*Tf zeros(lf,linf); zeros(linf,lf) eye(linf)]; 
B_r = [(Wf.')*(Bu-Z*Binf);(Winf.')*Binf];
C_r = [Cf*Tf (Cf*Y+Cu)*Tinf]; 


end

