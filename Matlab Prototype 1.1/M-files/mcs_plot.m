
%% this script is used to plot the reachable set of original model and reduced model for ISS Model
% It will get the data from file.gen in SpaceEx and then plot the data

% Author : Tran Hoang Dung
% Date : 21-5-2015
% Paper : Verification for large scalse system using order reduction method


% plot the reachable set of 5-order output abstraction of periodically switched synchronous motor control system 

% plot reachable set of y1-y2
fig = figure ;  plot_2d_vertices 'mcs_8_y1_y2.gen' 'b';
hold on;
V1 = [0.35 0.45; 0.35 0.6;0.4 0.45; 0.4 0.6];
P1 = polytope(V1);
plot(P1,'r');

