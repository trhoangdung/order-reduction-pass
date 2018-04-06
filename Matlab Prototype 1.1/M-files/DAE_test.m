E = [1 0 0 0; 0 0 1 0; 0 0 0 0; 0 0 0 0];
A = [0 1 0 0; 1 0 0 0; -1 0 0 1; 0 1 1 1];
B = [0; 0; 0; -1];
C= [0 0 1 0];
D = zeros(1,1);
sys = dss(A,B,C,D,E);

E1 = [0 1; 0 0];
A1 = eye(2); 
B1 = [-1;-1]; 
C1 = eye(2); 
sys1 = dss(A1,B1,C1,0,E1);

[m,n] = size(A); 
[re] = rank(E); 

% step 1: compute generalized Schur form
[EE,AA,VT,U] = qz(E,A,'real'); 
V = VT.'; 
Ef = EE(1:re,1:re);
Eu = EE(1:re,re+1:m);
Einf = EE(re+1:m,re+1:m);

BB = VT*B;
Bu = BB(1:re,:);
Binf = BB(re+1:m,:); 

CC = C*U;
Cf = CC(:,1:re);
Cu = CC(:,re+1:m);

Y = lyap((inv(Af))*Ef,-inv(Ainf)*Einf,-inv(Af)*(-Eu + Au*(inv(Ainf))*Einf));
Z = (Af*Y+Au)*(inv(Ainf));

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
%LinfT = chol(Xio);
%Linf = LinfT.'; 

r = rank(Lf.'*Ef*Rf);
%linf = rank(Linf.'*Ainf*Rinf);

[UU,ZZ,VV] = svd(Lf.'*Ef*Rf); 
%[U3,D3,V3] = svd(Linf.'*Ainf*Rinf);


