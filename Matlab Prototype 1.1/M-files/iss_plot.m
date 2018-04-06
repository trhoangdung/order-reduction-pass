plot_safespec(issr,'S','line')
plot_safespec(issr,'U','line')
plot_safespec(iss_sys,'S','line')
fig = figure(1);
plot_2d_vertices 'issr-10-y1-y2.gen' 'b';
fig = figure(2);
plot_2d_vertices 'issr-10-y1-y3.gen' 'b';
fig = figure(3);
plot_2d_vertices 'issr-10-y2-y3.gen' 'b';

fig = figure(4);
P = polytope(V); 
Pr = polytope(issr.S.matrix_A,-issr.S.matrix_B);
UPr = polytope(issr.U.matrix_A,-issr.U.matrix_B);
plot(P,'g',Pr,'b',UPr,'r');

