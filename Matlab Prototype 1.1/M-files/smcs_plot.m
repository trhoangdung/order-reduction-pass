
%% this script is used to plot the reachable set of original model and reduced model for ISS Model
% It will get the data from file.gen in SpaceEx and then plot the data

% Author : Tran Hoang Dung
% Date : 21-5-2015
% Paper : Verification for large scalse system using order reduction method


% plot the reachable set of 5-order output abstraction of periodically switched synchronous motor control system 

% plot reachable set of y1-y2
fig = figure ;  plot_2d_vertices 'smcs_yr1_yr2.gen' 'b';
hold on;

V1 = [0.25 0.12; 0.25 0.2;0.4 0.12; 0.4 0.2];
V2 = -V1; 

P1 = polytope(V1);
[E1,c1] = mpt_getInnerEllipsoid(P1); 
P2 = polytope(V2);
[E2,c2] = mpt_getInnerEllipsoid(P2); 

mpt_plotellip(E1,c1,'r');
mpt_plotellip(E2,c2,'r');

plot_safespec(mcsr_sub1,'U','line')
plot_safespec(mcsr_sub2,'U','line')




